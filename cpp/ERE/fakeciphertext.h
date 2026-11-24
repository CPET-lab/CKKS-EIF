#pragma once

#include <memory>
#include <utility>
#include <vector>
#include "error_bound.h"

namespace fct
{
    // Store one decoded nominal polynomial and one error-radius bound per input interval.
    class Ciphertext
    {
    public:
        std::vector<double> g;
        std::vector<double> r;
        std::shared_ptr<const std::vector<std::pair<double, double>>> intervals;
        double scale;
        int depth;

        // Create an empty state for an intermediate result.
        Ciphertext();
        // Initialize repeated singleton inputs using the legacy scalar constructor arguments.
        Ciphertext(double x, int size, double err = 0.0, double scale = 1.0, int depth = 0);
        // Partition an input domain and initialize its nominal input and error radii.
        Ciphertext(double a, double b, double interval_size = 1.0 / 16384.0, double err = 0.0, double scale = 1.0, bool encode = true, int depth = 0);
        // Initialize a state on an explicit list of input intervals.
        Ciphertext(const std::vector<std::pair<double, double>>& intervals, double err = 0.0, double scale = 1.0, bool encode = true, int depth = 0);
        // Interpret the supplied points as consecutive interval boundaries.
        Ciphertext(const std::vector<double>& x, double err = 0.0, double scale = 1.0, bool encode = true, int depth = 0);
        // Reset every radius using an error expressed in the legacy encoded units.
        void set_error(double err);
    };

    // Generate endpoint-inclusive interval boundaries.
    std::vector<double> generate_points(double a, double b, double interval_size);
    // Bound the departure of a polynomial from its endpoint interpolant on an interval.
    double kappa(const std::vector<double>& coeff, std::pair<double, double> interval);
    // Bound the absolute value of a polynomial throughout an interval.
    double max_abs(const std::vector<double>& coeff, std::pair<double, double> interval);

    // Add nominal polynomials and their interval-wise error radii.
    Ciphertext add(const Ciphertext& ct1, const Ciphertext& ct2, ErrBound& eb);
    // Add a plaintext constant and include its non-integer encoding error.
    Ciphertext add_pt(const Ciphertext& ct, double scalar);
    // Multiply nominal polynomials and propagate interval-wise multiplication errors.
    Ciphertext mult(const Ciphertext& ct1, const Ciphertext& ct2, ErrBound& eb);
    // Square a nominal polynomial and propagate interval-wise squaring errors.
    Ciphertext square(const Ciphertext& ct, ErrBound& eb);
    // Multiply by a scalar and include encoding and rescaling error for non-integer coefficients.
    Ciphertext mult_pt(double scalar, const Ciphertext& ct, ErrBound& eb, bool rescale = true);

    // Return decoded upper and lower bounds after curvature and error-radius correction.
    std::pair<std::vector<double>, std::vector<double>> dec(const Ciphertext& ct);
    // Aggregate all corrected interval bounds into a single output interval.
    std::pair<double, double> dec_interval(const Ciphertext& ct);
    // Return the corrected output interval for each input interval.
    std::vector<std::pair<double, double>> dec_intervals(const Ciphertext& ct);

    // Execute an existing EvalStep using nominal polynomials and interval error radii.
    Ciphertext evaluate(const Ciphertext& x, const EvalStep& es, ErrBound& eb);
    // Execute a composition without resetting its nominal polynomial or accumulated radii.
    Ciphertext evaluate_composite(const Ciphertext& x, const std::vector<EvalStep>& eval_steps, ErrBound& eb);
    // Estimate a polynomial's range with the existing call signature and zero initial error.
    std::pair<double, double> estimate_range(ErrBound& eb, std::pair<double, double> interval, const EvalStep& es, double interval_size);
    // Estimate a composition's range on the original partition with zero initial error.
    std::pair<double, double> estimate_composite_range(ErrBound& eb, std::pair<double, double> interval, const std::vector<EvalStep>& eval_steps, double interval_size);

}
