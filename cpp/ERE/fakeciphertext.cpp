#include "fakeciphertext.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fct
{
    namespace
    {
        // Reject non-finite or nonpositive scale metadata.
        void validate_scale(double scale)
        {
            if(!std::isfinite(scale) || scale <= 0.0)
                throw std::invalid_argument("fct: scale must be finite and positive");
        }

        // Order finite interval endpoints without changing the identity of the interval.
        std::pair<double, double> normalize_interval(std::pair<double, double> interval)
        {
            if(!std::isfinite(interval.first) || !std::isfinite(interval.second))
                throw std::invalid_argument("fct: interval endpoints must be finite");
            if(interval.first > interval.second)
                std::swap(interval.first, interval.second);
            return interval;
        }

        // Check that a nominal polynomial has finite coefficients and a defined constant term.
        void validate_polynomial(const std::vector<double>& coeff)
        {
            if(coeff.empty())
                throw std::invalid_argument("fct: empty nominal polynomial");
            for(double value: coeff)
                if(!std::isfinite(value))
                    throw std::invalid_argument("fct: polynomial coefficients must be finite");
        }

        // Convert a finite intermediate to double without rounding an upper bound inward.
        double upper_bound(long double value)
        {
            double result = static_cast<double>(value);
            if(!std::isfinite(value) || !std::isfinite(result))
                throw std::overflow_error("fct: non-finite bound");
            if(static_cast<long double>(result) < value)
                result = std::nextafter(result, std::numeric_limits<double>::infinity());
            return result;
        }

        // Convert a finite intermediate to double without rounding a lower bound inward.
        double lower_bound(long double value)
        {
            double result = static_cast<double>(value);
            if(!std::isfinite(value) || !std::isfinite(result))
                throw std::overflow_error("fct: non-finite bound");
            if(static_cast<long double>(result) > value)
                result = std::nextafter(result, -std::numeric_limits<double>::infinity());
            return result;
        }

        // Evaluate ascending-order polynomial coefficients using Horner's rule.
        long double polynomial_value(const std::vector<double>& coeff, long double x)
        {
            long double value = 0.0L;
            for(auto coefficient = coeff.rbegin(); coefficient != coeff.rend(); ++coefficient)
                value = value * x + *coefficient;
            return value;
        }

        // Remove only exactly zero leading coefficients while retaining a zero constant polynomial.
        std::vector<double> trim_coefficients(std::vector<double> coeff)
        {
            while(coeff.size() > 1 && coeff.back() == 0.0)
                coeff.pop_back();
            validate_polynomial(coeff);
            return coeff;
        }

        // Add ascending-order coefficient vectors with zero padding.
        std::vector<double> add_polynomials(const std::vector<double>& a, const std::vector<double>& b)
        {
            std::vector<double> result(std::max(a.size(), b.size()), 0.0);
            for(std::size_t i = 0; i < a.size(); ++i)
                result[i] += a[i];
            for(std::size_t i = 0; i < b.size(); ++i)
                result[i] += b[i];
            return trim_coefficients(std::move(result));
        }

        // Multiply ascending-order coefficient vectors by convolution.
        std::vector<double> multiply_polynomials(const std::vector<double>& a, const std::vector<double>& b)
        {
            std::vector<double> result(a.size() + b.size() - 1, 0.0);
            for(std::size_t i = 0; i < a.size(); ++i)
                for(std::size_t j = 0; j < b.size(); ++j)
                    result[i + j] += a[i] * b[j];
            return trim_coefficients(std::move(result));
        }

        // Compute the curvature correction from midpoint Taylor coefficients without factorials.
        long double curvature(const std::vector<double>& coeff, std::pair<double, double> interval)
        {
            if(coeff.size() <= 2 || interval.first == interval.second)
                return 0.0L;
            const long double p = interval.first;
            const long double q = interval.second;
            const long double midpoint = (p + q) / 2.0L;
            const long double radius = (q - p) / 2.0L;
            std::vector<long double> shifted(coeff.size(), 0.0L);
            for(auto coefficient = coeff.rbegin(); coefficient != coeff.rend(); ++coefficient)
            {
                for(std::size_t i = shifted.size() - 1; i > 0; --i)
                    shifted[i] = midpoint * shifted[i] + shifted[i - 1];
                shifted[0] = midpoint * shifted[0] + *coefficient;
            }
            long double derivative_bound = 0.0L;
            for(std::size_t i = shifted.size(); i-- > 2; )
                derivative_bound = derivative_bound * radius
                    + static_cast<long double>(i) * (i - 1) * std::abs(shifted[i]);
            return radius * radius * derivative_bound / 2.0L;
        }

        // Verify that each interval has a finite nonnegative radius and a common nominal polynomial.
        void validate_ciphertext(const Ciphertext& ct)
        {
            validate_scale(ct.scale);
            validate_polynomial(ct.g);
            if(!ct.intervals || ct.intervals->empty() || ct.r.size() != ct.intervals->size())
                throw std::invalid_argument("fct: interval and radius sizes must match and be nonempty");
            for(double radius: ct.r)
                if(!std::isfinite(radius) || radius < 0.0)
                    throw std::invalid_argument("fct: error radii must be finite and nonnegative");
        }

        // Require binary operands to use the same input partition.
        void validate_operands(const Ciphertext& ct1, const Ciphertext& ct2)
        {
            validate_ciphertext(ct1);
            validate_ciphertext(ct2);
            if(ct1.intervals != ct2.intervals && *ct1.intervals != *ct2.intervals)
                throw std::invalid_argument("fct: operands must share the input partition");
        }

        // Allocate an intermediate state while sharing the immutable input partition.
        Ciphertext result_state(const Ciphertext& input, std::vector<double> coeff, int depth)
        {
            Ciphertext result;
            result.g = std::move(coeff);
            result.r.resize(input.r.size(), 0.0);
            result.intervals = input.intervals;
            result.scale = input.scale;
            result.depth = depth;
            return result;
        }

        // Integer constants are exact under the adopted scalar-encoding model.
        long double encoding_error(double scalar, double scale)
        {
            return scalar == std::trunc(scalar) ? 0.0L : 0.5L / scale;
        }

        // Convert the existing encoded rescaling bound to decoded units.
        double rescaling_error(const ErrBound& eb)
        {
            validate_scale(eb.scale);
            if(!std::isfinite(eb.Bs) || eb.Bs < 0.0)
                throw std::invalid_argument("fct: rescaling error must be finite and nonnegative");
            return upper_bound(static_cast<long double>(eb.Bs) / eb.scale);
        }

        // Build consecutive input intervals from a finite sorted boundary vector.
        std::vector<std::pair<double, double>> intervals_from_points(std::vector<double> points)
        {
            if(points.empty())
                throw std::invalid_argument("fct: empty boundary vector");
            for(double point: points)
                if(!std::isfinite(point))
                    throw std::invalid_argument("fct: boundaries must be finite");
            std::sort(points.begin(), points.end());
            points.erase(std::unique(points.begin(), points.end()), points.end());
            if(points.size() == 1)
                return {{points.front(), points.front()}};
            std::vector<std::pair<double, double>> intervals;
            intervals.reserve(points.size() - 1);
            for(std::size_t i = 1; i < points.size(); ++i)
                intervals.emplace_back(points[i - 1], points[i]);
            return intervals;
        }

        // Initialize decoded input coordinates and convert the legacy encoded error units.
        void initialize(Ciphertext& ct, std::vector<std::pair<double, double>> intervals, double err, double scale, bool encode, int depth)
        {
            validate_scale(scale);
            if(intervals.empty() || depth < 0)
                throw std::invalid_argument("fct: nonempty intervals and nonnegative depth are required");
            for(auto& interval: intervals)
            {
                interval = normalize_interval(interval);
                if(!encode)
                {
                    interval.first /= scale;
                    interval.second /= scale;
                    interval = normalize_interval(interval);
                }
            }
            ct.g = {0.0, 1.0};
            ct.intervals = std::make_shared<const std::vector<std::pair<double, double>>>(std::move(intervals));
            ct.scale = scale;
            ct.depth = depth;
            ct.set_error(err);
        }
    }

    // Create an empty state for an intermediate result.
    Ciphertext::Ciphertext(): scale(1.0), depth(0)
    {
    }

    // Initialize repeated singleton inputs with the existing scalar constructor convention.
    Ciphertext::Ciphertext(double x, int size, double err, double scale, int depth)
    {
        if(size <= 0)
            throw std::invalid_argument("fct: size must be positive");
        initialize(*this, std::vector<std::pair<double, double>>(size, {x, x}), err, scale, true, depth);
    }

    // Partition an input domain into endpoint-inclusive intervals.
    Ciphertext::Ciphertext(double a, double b, double interval_size, double err, double scale, bool encode, int depth)
    {
        initialize(*this, intervals_from_points(generate_points(a, b, interval_size)), err, scale, encode, depth);
    }

    // Initialize the nominal input and radii on explicit input intervals.
    Ciphertext::Ciphertext(const std::vector<std::pair<double, double>>& intervals, double err, double scale, bool encode, int depth)
    {
        initialize(*this, intervals, err, scale, encode, depth);
    }

    // Use the supplied points as input-interval boundaries rather than isolated samples.
    Ciphertext::Ciphertext(const std::vector<double>& x, double err, double scale, bool encode, int depth)
    {
        initialize(*this, intervals_from_points(x), err, scale, encode, depth);
    }

    // Replace the accumulated radii with one decoded initial-error bound.
    void Ciphertext::set_error(double err)
    {
        validate_scale(scale);
        if(!std::isfinite(err) || err < 0.0 || !intervals)
            throw std::invalid_argument("fct: error must be finite and nonnegative on a valid partition");
        r.assign(intervals->size(), upper_bound(static_cast<long double>(err) / scale));
    }

    // Copy the endpoint-inclusive midpoint sampling algorithm to generate partition boundaries.
    std::vector<double> generate_points(double a, double b, double interval_size)
    {
        if(!std::isfinite(a) || !std::isfinite(b))
            throw std::invalid_argument("generate_points: interval endpoints must be finite");
        if(!std::isfinite(interval_size) || interval_size <= 0.0)
            throw std::invalid_argument("generate_points: interval_size must be positive");
        if(a > b)
            std::swap(a, b);
        if(a == b)
            return {a};
        if(!std::isfinite(b - a))
            throw std::invalid_argument("generate_points: interval width must be finite");

        const double midpoint = std::midpoint(a, b);
        std::vector<double> left_points;
        std::vector<double> right_points;
        double previous_left = midpoint;
        double previous_right = midpoint;
        for(std::size_t step = 1; ; ++step)
        {
            const double offset = interval_size * static_cast<double>(step);
            const double left = midpoint - offset;
            const double right = midpoint + offset;
            if(left <= a && right >= b)
                break;
            if(left > a)
            {
                if(left >= previous_left)
                    throw std::invalid_argument("generate_points: interval_size is below boundary precision");
                left_points.push_back(left);
                previous_left = left;
            }
            if(right < b)
            {
                if(right <= previous_right)
                    throw std::invalid_argument("generate_points: interval_size is below boundary precision");
                right_points.push_back(right);
                previous_right = right;
            }
        }
        std::vector<double> points;
        points.reserve(left_points.size() + right_points.size() + 3);
        points.push_back(a);
        points.insert(points.end(), left_points.rbegin(), left_points.rend());
        points.push_back(midpoint);
        points.insert(points.end(), right_points.begin(), right_points.end());
        points.push_back(b);
        return points;
    }

    // Compute the manuscript's curvature correction on a finite ordered interval.
    double kappa(const std::vector<double>& coeff, std::pair<double, double> interval)
    {
        validate_polynomial(coeff);
        interval = normalize_interval(interval);
        return upper_bound(curvature(coeff, interval));
    }

    // Bound a nominal polynomial's magnitude using endpoint values and curvature.
    double max_abs(const std::vector<double>& coeff, std::pair<double, double> interval)
    {
        validate_polynomial(coeff);
        interval = normalize_interval(interval);
        return upper_bound(std::max(std::abs(polynomial_value(coeff, interval.first)),
                                    std::abs(polynomial_value(coeff, interval.second)))
                           + curvature(coeff, interval));
    }

    // Add nominal coefficient vectors and accumulate the operands' interval radii.
    Ciphertext add(const Ciphertext& ct1, const Ciphertext& ct2, ErrBound&)
    {
        validate_operands(ct1, ct2);
        Ciphertext result = result_state(ct1, add_polynomials(ct1.g, ct2.g), std::max(ct1.depth, ct2.depth));
        for(std::size_t i = 0; i < result.r.size(); ++i)
            result.r[i] = upper_bound(static_cast<long double>(ct1.r[i]) + ct2.r[i]);
        return result;
    }

    // Include encoding error for a non-integer plaintext addend.
    Ciphertext add_pt(const Ciphertext& ct, double scalar)
    {
        validate_ciphertext(ct);
        if(!std::isfinite(scalar))
            throw std::invalid_argument("fct: scalar must be finite");
        Ciphertext result = ct;
        result.g[0] += scalar;
        result.g = trim_coefficients(std::move(result.g));
        const long double scalar_error = encoding_error(scalar, ct.scale);
        for(double& radius: result.r)
            radius = upper_bound(static_cast<long double>(radius) + scalar_error);
        return result;
    }

    // Multiply nominal polynomials and apply the interval magnitude bounds to the error formula.
    Ciphertext mult(const Ciphertext& ct1, const Ciphertext& ct2, ErrBound& eb)
    {
        validate_operands(ct1, ct2);
        const double error = rescaling_error(eb);
        Ciphertext result = result_state(ct1, multiply_polynomials(ct1.g, ct2.g), std::max(ct1.depth, ct2.depth) + 1);
        for(std::size_t i = 0; i < result.r.size(); ++i)
        {
            const long double m1 = max_abs(ct1.g, (*ct1.intervals)[i]);
            const long double m2 = max_abs(ct2.g, (*ct2.intervals)[i]);
            const long double r1 = ct1.r[i];
            const long double r2 = ct2.r[i];
            result.r[i] = upper_bound(m1 * r2 + m2 * r1 + r1 * r2 + error);
        }
        return result;
    }

    // Square the common nominal polynomial and propagate the corresponding interval error bound.
    Ciphertext square(const Ciphertext& ct, ErrBound& eb)
    {
        validate_ciphertext(ct);
        const double error = rescaling_error(eb);
        Ciphertext result = result_state(ct, multiply_polynomials(ct.g, ct.g), ct.depth + 1);
        for(std::size_t i = 0; i < result.r.size(); ++i)
        {
            const long double magnitude = max_abs(ct.g, (*ct.intervals)[i]);
            const long double radius = ct.r[i];
            result.r[i] = upper_bound(2.0L * magnitude * radius + radius * radius + error);
        }
        return result;
    }

    // Scale the nominal polynomial and propagate scalar encoding and rescaling errors per interval.
    Ciphertext mult_pt(double scalar, const Ciphertext& ct, ErrBound& eb, bool rescale)
    {
        validate_ciphertext(ct);
        if(!std::isfinite(scalar))
            throw std::invalid_argument("fct: scalar must be finite");
        const bool needs_rescale = rescale && scalar != std::trunc(scalar);
        const double error = needs_rescale ? rescaling_error(eb) : 0.0;
        const long double scalar_error = needs_rescale ? 0.5L / ct.scale : 0.0L;
        std::vector<double> coeff = ct.g;
        for(double& value: coeff)
            value *= scalar;
        Ciphertext result = result_state(ct, trim_coefficients(std::move(coeff)), ct.depth + (needs_rescale ? 1 : 0));
        for(std::size_t i = 0; i < result.r.size(); ++i)
        {
            const long double radius = ct.r[i];
            const long double magnitude = scalar_error == 0.0L
                ? 0.0L : max_abs(ct.g, (*ct.intervals)[i]);
            // The encoded scalar also multiplies the operand's existing error.
            result.r[i] = upper_bound(std::abs(static_cast<long double>(scalar)) * radius
                + scalar_error * (magnitude + radius) + error);
        }
        return result;
    }

    // Form decoded interval bounds by adding both nominal curvature and accumulated CKKS error.
    std::pair<std::vector<double>, std::vector<double>> dec(const Ciphertext& ct)
    {
        validate_ciphertext(ct);
        std::vector<double> high, low;
        high.reserve(ct.r.size());
        low.reserve(ct.r.size());
        for(std::size_t i = 0; i < ct.r.size(); ++i)
        {
            const auto interval = (*ct.intervals)[i];
            const long double left = polynomial_value(ct.g, interval.first);
            const long double right = polynomial_value(ct.g, interval.second);
            const long double correction = curvature(ct.g, interval) + ct.r[i];
            low.push_back(lower_bound(std::min(left, right) - correction));
            high.push_back(upper_bound(std::max(left, right) + correction));
        }
        return {high, low};
    }

    // Aggregate the corrected intervals into the final output range.
    std::pair<double, double> dec_interval(const Ciphertext& ct)
    {
        auto [high, low] = dec(ct);
        return {*std::min_element(low.begin(), low.end()), *std::max_element(high.begin(), high.end())};
    }

    // Return one corrected lower/upper range per input interval.
    std::vector<std::pair<double, double>> dec_intervals(const Ciphertext& ct)
    {
        auto [high, low] = dec(ct);
        std::vector<std::pair<double, double>> intervals;
        intervals.reserve(high.size());
        for(std::size_t i = 0; i < high.size(); ++i)
            intervals.emplace_back(low[i], high[i]);
        return intervals;
    }

    // Evaluate the sequence using shared powers and interval error radii.
    Ciphertext evaluate(const Ciphertext& x, const EvalStep& es, ErrBound& eb)
    {
        validate_ciphertext(x);
        if(es.eval_step.empty())
            throw std::invalid_argument("evaluate: empty evaluation sequence");
        std::map<std::string, Ciphertext> fpowers;
        std::map<std::string, double> coeffs;
        std::map<std::string, Ciphertext> fterms;
        fpowers.insert_or_assign("P1", x);
        Ciphertext res = x;

        // Resolve a stored power or intermediate using the existing key prefixes.
        const auto operand = [&](const std::string& key) -> const Ciphertext& {
            if(key.empty())
                throw std::invalid_argument("evaluate: empty operand key");
            if(key[0] == 'P')
                return fpowers.at(key);
            if(key[0] == 'T')
                return fterms.at(key);
            throw std::invalid_argument("evaluate: invalid ciphertext key " + key);
        };

        for(const auto& step: es.eval_step)
        {
            if(step.key1.empty())
                throw std::invalid_argument("evaluate: empty first operand key");
            switch(step.op)
            {
                case 'o':
                    if(step.key1[0] == 'P')
                    {
                        res = operand(step.key1);
                        continue;
                    }
                    if(step.key1[0] != 'C')
                        throw std::invalid_argument("evaluate: invalid constant key " + step.key1);
                    coeffs.insert_or_assign(step.key1, std::stod(step.key2));
                    res = result_state(x, {coeffs.at(step.key1)}, 0);
                    validate_polynomial(res.g);
                    std::fill(res.r.begin(), res.r.end(), upper_bound(encoding_error(res.g[0], res.scale)));
                    break;

                case '+':
                case 'x':
                {
                    if(step.key2.empty() || step.save_key.empty())
                        throw std::invalid_argument("evaluate: empty operand or result key");
                    const bool first_constant = step.key1[0] == 'C';
                    const bool second_constant = step.key2[0] == 'C';
                    if(first_constant && second_constant)
                    {
                        const double a = coeffs.at(step.key1);
                        const double b = coeffs.at(step.key2);
                        res = result_state(x, {step.op == '+' ? a + b : a * b}, 0);
                        validate_polynomial(res.g);
                        std::fill(res.r.begin(), res.r.end(), upper_bound(encoding_error(res.g[0], res.scale)));
                    }
                    else if(first_constant || second_constant)
                    {
                        const double scalar = coeffs.at(first_constant ? step.key1 : step.key2);
                        const Ciphertext& value = operand(first_constant ? step.key2 : step.key1);
                        res = step.op == '+' ? add_pt(value, scalar) : mult_pt(scalar, value, eb);
                    }
                    else if(step.op == '+')
                        res = add(operand(step.key1), operand(step.key2), eb);
                    else if(step.key1 == step.key2)
                        res = square(operand(step.key1), eb);
                    else
                        res = mult(operand(step.key1), operand(step.key2), eb);

                    switch(step.save_key[0])
                    {
                        case 'P':
                            fpowers.insert_or_assign(step.save_key, res);
                            break;
                        case 'T':
                            fterms.insert_or_assign(step.save_key, res);
                            break;
                        default:
                            throw std::invalid_argument("evaluate: invalid result key " + step.save_key);
                    }
                    break;
                }
                default:
                    throw std::invalid_argument("evaluate: unknown operation");
            }
        }
        return res;
    }

    // Carry the original input partition through every composition stage.
    Ciphertext evaluate_composite(const Ciphertext& x, const std::vector<EvalStep>& eval_steps, ErrBound& eb)
    {
        if(eval_steps.empty())
            throw std::invalid_argument("evaluate_composite: empty evaluation step list");
        Ciphertext res = x;
        for(const auto& es: eval_steps)
            res = evaluate(res, es, eb);
        return res;
    }

    // Preserve the existing range-estimation arguments while initializing every radius to zero.
    std::pair<double, double> estimate_range(ErrBound& eb, std::pair<double, double> interval, const EvalStep& es, double interval_size)
    {
        Ciphertext ct(interval.first, interval.second, interval_size, 0.0, eb.scale, true);
        ct = evaluate(ct, es, eb);
        return dec_interval(ct);
    }

    // Preserve the existing composite-range arguments and apply the final enclosure only once.
    std::pair<double, double> estimate_composite_range(ErrBound& eb, std::pair<double, double> interval, const std::vector<EvalStep>& eval_steps, double interval_size)
    {
        if(eval_steps.empty())
            throw std::invalid_argument("estimate_composite_range: empty evaluation step list");
        Ciphertext ct(interval.first, interval.second, interval_size, 0.0, eb.scale, true);
        ct = evaluate_composite(ct, eval_steps, eb);
        return dec_interval(ct);
    }
}
