#include "CKKS_params.h"
#include "error_bound.h"
#include "evaluate.h"
#include "fakeciphertext.h"
#include "json.hpp"
#include "multicore.h"
#include "parse.h"
#include "polyEval_class.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <syncstream>
#include <utility>
#include <vector>

namespace fs = std::filesystem;
using json = nlohmann::json;

namespace {

constexpr double kSigma = 3.2;
constexpr int kLogN = 17;
constexpr int kHammingWeight = 128;
constexpr int kScaleBits = 50;
constexpr int kPrecisionBits = 30;
constexpr double kExtremaScanSpacing = 1.0 / 65536.0;

struct Options {
    std::string domain_mode = "combined";
    int p_exponent = 2;
    int rounds = 10;
    int workers = 1;
    int ct_spacing_exponent = 16;
    int first_fct_exponent = 8;
    int last_fct_exponent = 12;
    fs::path output;
    fs::path coefficient_directory;
    fs::path ct_input;
};

struct Samples {
    std::vector<std::vector<double>> batches;
    std::size_t unique_count = 0;
};

struct PlaintextExtrema {
    double minimum_x = 0.0;
    double minimum_value = std::numeric_limits<double>::infinity();
    double maximum_x = 0.0;
    double maximum_value = -std::numeric_limits<double>::infinity();
    std::size_t critical_point_count = 0;
};

fs::path source_root()
{
#ifdef APPROX_SOURCE_DIR
    return fs::path(APPROX_SOURCE_DIR);
#else
    return fs::current_path();
#endif
}

void print_usage(const char* executable)
{
    std::cout
        << "Usage: " << executable << " [options]\n"
        << "  --p-exp N        Evaluate the p_num=N coefficient pipeline (default: 2)\n"
        << "  --rounds N       Number of independent CKKS encryptions (default: 10)\n"
        << "  --workers N      Number of parallel CT round workers (default: 1)\n"
        << "  --ct-exp N       CT reference spacing is 2^-N (default: 16)\n"
        << "  --fct-first N    First FCT spacing exponent (default: 8)\n"
        << "  --fct-last N     Last FCT spacing exponent (default: 12)\n"
        << "  --domain-mode M combined (default) or split input regions\n"
        << "  --output PATH    Output JSON path\n"
        << "  --coeff-dir PATH Coefficient directory (default: cache/ops/)\n"
        << "  --ct-input PATH Reuse CT bounds from JSON; only recompute Baseline and FCT\n"
        << "  --quick          One round with CT spacing 2^-12\n"
        << "  --help           Show this message\n";
}

int parse_positive_int(const std::string& value, const std::string& option)
{
    const int parsed = std::stoi(value);
    if(parsed <= 0)
        throw std::invalid_argument(option + " must be positive");
    return parsed;
}

Options parse_options(int argc, char** argv)
{
    Options options;

    for(int i = 1; i < argc; ++i)
    {
        const std::string argument = argv[i];
        auto next_value = [&](const std::string& option) -> std::string {
            if(++i >= argc)
                throw std::invalid_argument(option + " requires a value");
            return argv[i];
        };

        if(argument == "--domain-mode")
            options.domain_mode = next_value(argument);
        else if(argument == "--p-exp")
            options.p_exponent = parse_positive_int(next_value(argument), argument);
        else if(argument == "--rounds")
            options.rounds = parse_positive_int(next_value(argument), argument);
        else if(argument == "--workers" || argument == "--cores")
            options.workers = parse_positive_int(next_value(argument), argument);
        else if(argument == "--ct-exp")
            options.ct_spacing_exponent = parse_positive_int(next_value(argument), argument);
        else if(argument == "--fct-first")
            options.first_fct_exponent = parse_positive_int(next_value(argument), argument);
        else if(argument == "--fct-last")
            options.last_fct_exponent = parse_positive_int(next_value(argument), argument);
        else if(argument == "--output")
            options.output = fs::absolute(next_value(argument));
        else if(argument == "--coeff-dir")
            options.coefficient_directory = fs::absolute(next_value(argument));
        else if(argument == "--ct-input")
            options.ct_input = fs::absolute(next_value(argument));
        else if(argument == "--quick")
        {
            options.rounds = 1;
            options.ct_spacing_exponent = 12;
        }
        else if(argument == "--help")
        {
            print_usage(argv[0]);
            std::exit(0);
        }
        else
            throw std::invalid_argument("Unknown option: " + argument);
    }

    if(options.domain_mode != "combined" && options.domain_mode != "split")
        throw std::invalid_argument("--domain-mode must be combined or split");
    if(options.first_fct_exponent > options.last_fct_exponent)
        throw std::invalid_argument("--fct-first must not exceed --fct-last");
    if(options.p_exponent >= 31)
        throw std::invalid_argument("--p-exp must be less than 31");
    if(options.output.empty())
    {
        const std::string suffix = options.p_exponent == 2
            ? ""
            : "_p" + std::to_string(options.p_exponent);
        const std::string filename = std::string("ere_results")
            + (options.domain_mode == "combined" ? "_combined" : "") + suffix
            + ".json";
        options.output = source_root() / "test/ERE/result" / filename;
    }
    if(options.coefficient_directory.empty())
        options.coefficient_directory = source_root() / "cache/ops";
    if(!options.ct_input.empty() && options.ct_input == options.output)
        throw std::invalid_argument("CT input and output must be different files");

    return options;
}

std::vector<std::vector<double>> load_pipeline_coefficients(
    const Options& options)
{
    const std::string filename = remez_coeff_filename(
        options.p_exponent, kPrecisionBits, "depth"
    );
    std::vector<std::vector<double>> coefficients = parse_remez_coeff(
        options.p_exponent,
        kPrecisionBits,
        false,
        options.coefficient_directory.string(),
        filename
    );
    if(coefficients.empty())
        throw std::runtime_error(
            "No p_num=" + std::to_string(options.p_exponent)
            + ", e_num=" + std::to_string(kPrecisionBits)
            + " polynomials found in "
            + (options.coefficient_directory / filename).string()
        );
    return coefficients;
}

std::shared_ptr<Decomp> restore_decomposition(
    const std::vector<double>& coefficients,
    const json& cache)
{
    const std::string key = get_poly_type_key(coefficients);
    if(!cache.contains(key))
        throw std::runtime_error("Decomposition cache key not found: " + key);

    auto decomposition = reconstruct_decomp_from_cache(
        cache.at(key),
        Poly(coefficients),
        0
    );
    if(!decomposition)
        throw std::runtime_error("Failed to restore decomposition for key: " + key);
    return decomposition;
}

Samples build_samples(
    const std::vector<std::pair<double, double>>& domains,
    double spacing,
    std::size_t slot_count)
{
    Samples samples;
    for(const auto& domain : domains)
    {
        std::vector<double> points = fct::generate_points(
            domain.first,
            domain.second,
            spacing
        );
        samples.unique_count += points.size();

        for(std::size_t begin = 0; begin < points.size(); begin += slot_count)
        {
            const std::size_t end = std::min(begin + slot_count, points.size());
            std::vector<double> batch(points.begin() + begin, points.begin() + end);
            batch.resize(slot_count, domain.second);
            samples.batches.push_back(std::move(batch));
        }
    }
    return samples;
}

double evaluate_polynomial(
    const std::vector<double>& coefficients,
    double x)
{
    double value = 0.0;
    for(auto coefficient = coefficients.rbegin();
        coefficient != coefficients.rend();
        ++coefficient)
    {
        value = value * x + *coefficient;
    }
    return value;
}

double evaluate_polynomial_derivative(
    const std::vector<double>& coefficients,
    double x)
{
    double value = 0.0;
    for(std::size_t index = coefficients.size(); index-- > 1; )
        value = value * x + static_cast<double>(index) * coefficients[index];
    return value;
}

double evaluate_pipeline(
    const std::vector<std::vector<double>>& coefficients,
    double x)
{
    double value = x;
    for(const auto& polynomial : coefficients)
        value = evaluate_polynomial(polynomial, value);
    return value;
}

double evaluate_pipeline_derivative(
    const std::vector<std::vector<double>>& coefficients,
    double x)
{
    double value = x;
    double derivative = 1.0;
    for(const auto& polynomial : coefficients)
    {
        derivative *= evaluate_polynomial_derivative(polynomial, value);
        value = evaluate_polynomial(polynomial, value);
    }
    return derivative;
}

double bisect_critical_point(
    const std::vector<std::vector<double>>& coefficients,
    double left,
    double right)
{
    double left_value = evaluate_pipeline_derivative(coefficients, left);
    for(int iteration = 0; iteration < 100; ++iteration)
    {
        const double midpoint = left + (right - left) / 2.0;
        if(midpoint == left || midpoint == right)
            return midpoint;

        const double midpoint_value = evaluate_pipeline_derivative(
            coefficients,
            midpoint
        );
        if(midpoint_value == 0.0)
            return midpoint;

        if(std::signbit(left_value) == std::signbit(midpoint_value))
        {
            left = midpoint;
            left_value = midpoint_value;
        }
        else
            right = midpoint;
    }
    return left + (right - left) / 2.0;
}

PlaintextExtrema find_plaintext_extrema(
    const std::vector<std::pair<double, double>>& domains,
    const std::vector<std::vector<double>>& coefficients)
{
    PlaintextExtrema extrema;
    std::vector<double> critical_points;

    auto update_extrema = [&](double x) {
        const double value = evaluate_pipeline(coefficients, x);
        if(value < extrema.minimum_value)
        {
            extrema.minimum_x = x;
            extrema.minimum_value = value;
        }
        if(value > extrema.maximum_value)
        {
            extrema.maximum_x = x;
            extrema.maximum_value = value;
        }
    };

    auto add_critical_point = [&](double x) {
        for(const double existing : critical_points)
        {
            if(std::abs(existing - x) <= 1e-12)
                return;
        }
        critical_points.push_back(x);
        update_extrema(x);
    };

    for(const auto& domain : domains)
    {
        const double start = domain.first;
        const double end = domain.second;
        update_extrema(start);
        update_extrema(end);

        const std::size_t segment_count = std::max<std::size_t>(
            1,
            static_cast<std::size_t>(
                std::ceil((end - start) / kExtremaScanSpacing)
            )
        );
        double previous_x = start;
        double previous_derivative = evaluate_pipeline_derivative(
            coefficients,
            previous_x
        );

        for(std::size_t segment = 1; segment <= segment_count; ++segment)
        {
            const double x = start
                + (end - start)
                * static_cast<double>(segment)
                / static_cast<double>(segment_count);
            const double derivative = evaluate_pipeline_derivative(
                coefficients,
                x
            );

            if(previous_derivative == 0.0
                && previous_x > start
                && previous_x < end)
            {
                add_critical_point(previous_x);
            }
            if(derivative == 0.0 && x > start && x < end)
                add_critical_point(x);
            else if(std::isfinite(previous_derivative)
                && std::isfinite(derivative)
                && std::signbit(previous_derivative) != std::signbit(derivative))
            {
                add_critical_point(
                    bisect_critical_point(coefficients, previous_x, x)
                );
            }

            previous_x = x;
            previous_derivative = derivative;
        }
    }
    extrema.critical_point_count = critical_points.size();
    return extrema;
}

std::unique_ptr<CKKS_params> make_ckks_parameters(
    int required_depth,
    int p_exponent)
{
    std::vector<int> modulus_bits = {60};
    modulus_bits.insert(modulus_bits.end(), required_depth, kScaleBits);
    modulus_bits.push_back(60);

    return std::make_unique<CKKS_params>(
        modulus_bits,
        p_exponent,
        kPrecisionBits,
        kScaleBits,
        kLogN,
        kHammingWeight
    );
}

std::vector<std::pair<double, double>> evaluate_ct_ranges(
    const std::vector<Samples>& domain_samples,
    const std::vector<EvalStep>& evaluations,
    int rounds,
    int required_depth,
    int p_exponent,
    int workers)
{
    multicore::TaskExecutor executor(static_cast<std::size_t>(workers));
    std::vector<std::unique_ptr<CKKS_params>> worker_parameters(
        executor.worker_count()
    );

    const std::vector<std::vector<std::pair<double, double>>> round_ranges = executor.map(
        static_cast<std::size_t>(rounds),
        [&](std::size_t round, std::size_t worker_index) {
            std::unique_ptr<CKKS_params>& parameters = worker_parameters[worker_index];
            if(!parameters)
                parameters = make_ckks_parameters(required_depth, p_exponent);

            std::osyncstream(std::cout)
                << "    CT round " << round + 1 << '/' << rounds
                << " (worker " << worker_index + 1 << ")" << std::endl;

            std::vector<std::pair<double, double>> ranges(
                domain_samples.size(),
                {
                    std::numeric_limits<double>::infinity(),
                    -std::numeric_limits<double>::infinity(),
                }
            );
            for(std::size_t domain_index = 0;
                domain_index < domain_samples.size();
                ++domain_index)
            {
                auto& [lower, upper] = ranges[domain_index];
                for(const auto& plaintext : domain_samples[domain_index].batches)
                {
                    seal::Ciphertext encrypted = parameters->encrypt(plaintext);
                    seal::Ciphertext evaluated = encrypted;
                    for(const auto& evaluation : evaluations)
                    {
                        evaluated = ct::evaluate(
                            evaluated,
                            evaluation,
                            *parameters
                        );
                    }
                    std::vector<double> decoded = parameters->decode_ctxt(evaluated);

                    lower = std::min(
                        lower,
                        *std::min_element(decoded.begin(), decoded.end())
                    );
                    upper = std::max(
                        upper,
                        *std::max_element(decoded.begin(), decoded.end())
                    );
                }
            }
            return ranges;
        }
    );

    std::vector<std::pair<double, double>> ranges(
        domain_samples.size(),
        {
            std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
        }
    );
    for(const auto& round : round_ranges)
    {
        for(std::size_t domain_index = 0;
            domain_index < round.size();
            ++domain_index)
        {
            ranges[domain_index].first = std::min(
                ranges[domain_index].first,
                round[domain_index].first
            );
            ranges[domain_index].second = std::max(
                ranges[domain_index].second,
                round[domain_index].second
            );
        }
    }
    return ranges;
}

std::vector<std::pair<double, double>> evaluate_fct_ranges(
    ErrBound& error_bound,
    const std::vector<std::pair<double, double>>& domains,
    const std::vector<EvalStep>& evaluations,
    double spacing)
{
    if(evaluations.empty())
        throw std::invalid_argument("evaluate_fct_ranges: empty evaluation step list");
    std::vector<std::pair<double, double>> ranges;
    ranges.reserve(domains.size());
    for(const auto& domain : domains)
    {
        auto range = domain;
        // Match margin_test: re-encode each stage's output range with initial error zero.
        for(const auto& evaluation : evaluations)
        {
            range = fct::estimate_range(error_bound, range, evaluation, spacing);
        }
        ranges.push_back(range);
    }
    return ranges;
}

json range_json(const std::pair<double, double>& range)
{
    return json::array({range.first, range.second});
}

std::pair<double, double> combine_ranges(
    const std::vector<std::pair<double, double>>& ranges)
{
    double lower = std::numeric_limits<double>::infinity();
    double upper = -std::numeric_limits<double>::infinity();
    for(const auto& range : ranges)
    {
        lower = std::min(lower, range.first);
        upper = std::max(upper, range.second);
    }
    return {lower, upper};
}

PlaintextExtrema combine_extrema(
    const std::vector<PlaintextExtrema>& extrema_by_domain)
{
    PlaintextExtrema combined;
    for(const auto& extrema : extrema_by_domain)
    {
        if(extrema.minimum_value < combined.minimum_value)
        {
            combined.minimum_x = extrema.minimum_x;
            combined.minimum_value = extrema.minimum_value;
        }
        if(extrema.maximum_value > combined.maximum_value)
        {
            combined.maximum_x = extrema.maximum_x;
            combined.maximum_value = extrema.maximum_value;
        }
        combined.critical_point_count += extrema.critical_point_count;
    }
    return combined;
}

json extrema_json(const PlaintextExtrema& extrema)
{
    return {
        {"minimum", json::array({extrema.minimum_x, extrema.minimum_value})},
        {"maximum", json::array({extrema.maximum_x, extrema.maximum_value})},
        {"critical_point_count", extrema.critical_point_count},
        {"scan_spacing", kExtremaScanSpacing},
    };
}

json named_ranges_json(
    const std::vector<std::pair<double, double>>& ranges)
{
    if(ranges.size() == 1)
        return {{"combined", range_json(ranges[0])}};
    if(ranges.size() != 2)
        throw std::invalid_argument("Expected target and non-target ranges");
    return {
        {"target", range_json(ranges[0])},
        {"non_target", range_json(ranges[1])},
    };
}

json load_ct_result(const fs::path& path, const json& parameters, const json& stages)
{
    std::ifstream input(path);
    if(!input)
        throw std::runtime_error("Failed to open CT input file: " + path.string());
    json result;
    input >> result;
    const auto& saved_parameters = result.at("parameters");
    for(const char* key : {"p_exponent", "p", "precision_bits", "scale_bits",
        "log_n", "hamming_weight", "sigma", "clean_error", "rescale_error"})
    {
        if(saved_parameters.at(key) != parameters.at(key))
            throw std::runtime_error("CT input parameter mismatch: " + std::string(key));
    }
    const auto& saved_case = result.at("case");
    if(saved_case.at("stages") != stages)
        throw std::runtime_error("CT input polynomial stages or decomposition do not match");
    const auto saved_domains = saved_parameters.value(
        "ct_reference_domains", saved_parameters.at("evaluation_domains")
    );
    if(saved_domains.empty() || saved_case.at("ct_ranges").size() != saved_domains.size())
        throw std::runtime_error("CT reference domains and ranges are inconsistent");
    const auto& domains = parameters.at("evaluation_domains");
    if(saved_domains != domains)
    {
        if(domains.size() != 1 || saved_domains.empty())
            throw std::runtime_error("CT input evaluation domains do not match");
        for(const auto& domain : saved_domains)
        {
            if(!domain.is_array() || domain.size() != 2
                || domain.at(0).get<double>() < domains.at(0).at(0).get<double>()
                || domain.at(1).get<double>() > domains.at(0).at(1).get<double>()
                || domain.at(0).get<double>() > domain.at(1).get<double>())
                throw std::runtime_error("CT reference domain is outside the evaluation domain");
        }
    }
    std::vector<std::pair<double, double>> ranges;
    for(const auto& [domain, range] : saved_case.at("ct_ranges").items())
    {
        if(!range.is_array() || range.size() != 2)
            throw std::runtime_error("Invalid CT input range: " + domain);
        const double low = range.at(0).get<double>();
        const double high = range.at(1).get<double>();
        if(!std::isfinite(low) || !std::isfinite(high) || low > high)
            throw std::runtime_error("Invalid CT input bounds: " + domain);
        ranges.emplace_back(low, high);
    }
    if(saved_case.at("ct_range") != range_json(combine_ranges(ranges)))
        throw std::runtime_error("CT input combined range does not match domain ranges");
    return result;
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const Options options = parse_options(argc, argv);
        const fs::path root = source_root();

        ErrBound error_bound(kSigma, 1 << kLogN, kHammingWeight, kScaleBits);
        const double clean_error = error_bound.Bc / error_bound.scale;
        const int p = 1 << options.p_exponent;
        const std::vector<std::pair<double, double>> evaluation_domains =
            options.domain_mode == "combined"
            ? std::vector<std::pair<double, double>>{{0.0, static_cast<double>(p - 1) + clean_error}}
            : std::vector<std::pair<double, double>>{
                {0.0, clean_error},
                {1.0 - clean_error, static_cast<double>(p - 1) + clean_error},
            };
        json domains_json = json::array();
        for(const auto& domain : evaluation_domains)
            domains_json.push_back(range_json(domain));

        const fs::path cache_path = root / "cache/opd/decomp_cache.json";
        const json decomposition_cache = load_decomp_cache(cache_path.string());
        if(decomposition_cache.is_null())
            throw std::runtime_error("Failed to load decomposition cache: " + cache_path.string());

        json result;
        result["parameters"] = {
            {"coefficient_file", (options.coefficient_directory / remez_coeff_filename(
                options.p_exponent, kPrecisionBits, "depth"
            )).string()},
            {"p_exponent", options.p_exponent},
            {"p", p},
            {"precision_bits", kPrecisionBits},
            {"scale_bits", kScaleBits},
            {"log_n", kLogN},
            {"hamming_weight", kHammingWeight},
            {"sigma", kSigma},
            {"rounds", options.rounds},
            {"workers_requested", options.workers},
            {"workers_used", std::min(options.workers, options.rounds)},
            {"ct_spacing_exponent", options.ct_spacing_exponent},
            {"ct_spacing", std::ldexp(1.0, -options.ct_spacing_exponent)},
            {"fct_first_exponent", options.first_fct_exponent},
            {"fct_last_exponent", options.last_fct_exponent},
            {"clean_error", clean_error},
            {"rescale_error", error_bound.Bs / error_bound.scale},
            {"domain_mode", options.domain_mode},
            {"evaluation_domains", domains_json},
            {"ct_reference_domains", domains_json},
            {"ct_reference_status", "fresh_same_domain"},
        };

        const std::vector<std::vector<double>> pipeline_coefficients =
            load_pipeline_coefficients(options);
        std::vector<EvalStep> evaluations;
        json stages = json::array();
        int total_depth = 0;
        int composition_degree = 1;

        for(std::size_t index = 0; index < pipeline_coefficients.size(); ++index)
        {
            const std::vector<double>& coefficients = pipeline_coefficients[index];
            const std::shared_ptr<Decomp> decomposition = restore_decomposition(
                coefficients,
                decomposition_cache
            );
            evaluations.emplace_back(decomposition);
            const int degree = static_cast<int>(coefficients.size()) - 1;
            total_depth += decomposition->comp.depth;
            composition_degree *= degree;
            stages.push_back({
                {"index", index + 1},
                {"degree", degree},
                {"coefficients", coefficients},
                {"decomposition_depth", decomposition->comp.depth},
            });
        }

        json ct_result;
        if(!options.ct_input.empty())
        {
            ct_result = load_ct_result(options.ct_input, result["parameters"], stages);
            for(const char* key : {"rounds", "workers_requested", "workers_used",
                "ct_spacing_exponent", "ct_spacing"})
                result["parameters"][key] = ct_result.at("parameters").at(key);
            result["parameters"]["ct_input"] = options.ct_input.string();
            const auto& saved_parameters = ct_result.at("parameters");
            const auto reference_domains = saved_parameters.value(
                "ct_reference_domains", saved_parameters.at("evaluation_domains")
            );
            result["parameters"]["ct_reference_domains"] = reference_domains;
            result["parameters"]["ct_reference_status"] = reference_domains == domains_json
                ? "reused_same_domain"
                : "reused_original_subset; no_CT_on_added_inputs";
        }

        const double ct_spacing = std::ldexp(1.0, -options.ct_spacing_exponent);
        std::vector<Samples> ct_samples_by_domain;
        std::vector<PlaintextExtrema> plaintext_extrema_by_domain;
        ct_samples_by_domain.reserve(evaluation_domains.size());
        plaintext_extrema_by_domain.reserve(evaluation_domains.size());
        for(const auto& domain : evaluation_domains)
        {
            if(options.ct_input.empty())
                ct_samples_by_domain.push_back(build_samples(
                    {domain},
                    ct_spacing,
                    std::size_t{1} << (kLogN - 1)
                ));
            plaintext_extrema_by_domain.push_back(find_plaintext_extrema(
                {domain},
                pipeline_coefficients
            ));
        }
        const PlaintextExtrema plaintext_extrema = combine_extrema(
            plaintext_extrema_by_domain
        );
        const double rescale_error = error_bound.Bs / error_bound.scale;
        std::vector<std::pair<double, double>> baseline_ranges;
        baseline_ranges.reserve(plaintext_extrema_by_domain.size());
        for(const auto& extrema : plaintext_extrema_by_domain)
        {
            baseline_ranges.emplace_back(
                extrema.minimum_value - rescale_error,
                extrema.maximum_value + rescale_error
            );
        }
        const std::pair<double, double> baseline_range = combine_ranges(
            baseline_ranges
        );

        std::size_t ct_unique_samples = 0;
        std::size_t ct_batches_per_round = 0;
        for(const auto& samples : ct_samples_by_domain)
        {
            ct_unique_samples += samples.unique_count;
            ct_batches_per_round += samples.batches.size();
        }

        std::cout << "ERE, composite degrees ";
        for(std::size_t index = 0; index < stages.size(); ++index)
        {
            if(index > 0)
                std::cout << " -> ";
            std::cout << stages[index]["degree"].get<int>();
        }
        std::cout << ", total depth " << total_depth << std::endl;
        std::vector<std::pair<double, double>> ct_ranges;
        if(options.ct_input.empty())
            ct_ranges = evaluate_ct_ranges(
                ct_samples_by_domain,
                evaluations,
                options.rounds,
                std::max(1, total_depth),
                options.p_exponent,
                options.workers
            );
        else
        {
            std::cout << "    Reusing CT bounds from " << options.ct_input << std::endl;
            const std::vector<std::string> domain_names =
                ct_result.at("case").at("ct_ranges").contains("combined")
                ? std::vector<std::string>{"combined"}
                : std::vector<std::string>{"target", "non_target"};
            for(const auto& domain : domain_names)
            {
                const auto& range = ct_result.at("case").at("ct_ranges").at(domain);
                ct_ranges.emplace_back(range.at(0).get<double>(), range.at(1).get<double>());
            }
        }
        const auto ct_range = combine_ranges(ct_ranges);

        json case_result = {
            {"composition_degree", composition_degree},
            {"total_decomposition_depth", total_depth},
            {"stages", stages},
            {"ct_unique_samples", ct_unique_samples},
            {"ct_batches_per_round", ct_batches_per_round},
            {"plaintext_extrema", extrema_json(plaintext_extrema)},
            {"baseline_method", "post_composition_rescale_expansion"},
            {"baseline_range", range_json(baseline_range)},
            {"baseline_ranges", named_ranges_json(baseline_ranges)},
            {"ct_range", range_json(ct_range)},
            {"ct_ranges", named_ranges_json(ct_ranges)},
            {"fct_method", "stagewise_range_reencoding"},
            {"fct_ranges", json::array()},
        };
        if(evaluation_domains.size() == 1)
            case_result["plaintext_extrema_by_domain"] = {
                {"combined", extrema_json(plaintext_extrema_by_domain[0])},
            };
        else
            case_result["plaintext_extrema_by_domain"] = {
                {"target", extrema_json(plaintext_extrema_by_domain[0])},
                {"non_target", extrema_json(plaintext_extrema_by_domain[1])},
            };
        if(options.ct_input.empty())
        {
            if(evaluation_domains.size() == 1)
                case_result["ct_unique_samples_by_domain"] = {
                    {"combined", ct_samples_by_domain[0].unique_count},
                };
            else
                case_result["ct_unique_samples_by_domain"] = {
                    {"target", ct_samples_by_domain[0].unique_count},
                    {"non_target", ct_samples_by_domain[1].unique_count},
                };
        }
        else
        {
            for(const char* key : {"ct_range", "ct_ranges", "ct_unique_samples",
                "ct_batches_per_round", "ct_unique_samples_by_domain"})
                case_result[key] = ct_result.at("case").at(key);
        }

        for(int exponent = options.first_fct_exponent;
            exponent <= options.last_fct_exponent;
            ++exponent)
        {
            const double spacing = std::ldexp(1.0, -exponent);
            const auto fct_ranges = evaluate_fct_ranges(
                error_bound,
                evaluation_domains,
                evaluations,
                spacing
            );
            const auto fct_range = combine_ranges(fct_ranges);
            std::cout << "    FCT h=2^-" << exponent;
            for(std::size_t index = 0; index < fct_ranges.size(); ++index)
            {
                const auto reference = fct_ranges.size() == 1 ? ct_range : ct_ranges[index];
                std::cout << (fct_ranges.size() == 1 ? " combined" : (index == 0 ? " target" : " non_target"))
                    << " lower=" << (fct_ranges[index].first <= reference.first ? "PASS" : "FAIL")
                    << " upper=" << (fct_ranges[index].second >= reference.second ? "PASS" : "FAIL");
            }
            std::cout << std::endl;
            case_result["fct_ranges"].push_back({
                {"exponent", exponent},
                {"spacing", spacing},
                {"range", range_json(fct_range)},
                {"ranges", named_ranges_json(fct_ranges)},
            });
        }
        result["case"] = std::move(case_result);

        fs::create_directories(options.output.parent_path());
        std::ofstream output(options.output);
        if(!output)
            throw std::runtime_error("Failed to open output file: " + options.output.string());
        output << std::setprecision(17) << result.dump(2) << '\n';
        output.close();
        std::cout << "Wrote " << options.output << std::endl;

        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << "ERE experiment failed: " << error.what() << std::endl;
    }
    return 1;
}
