#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "parse.h"
#include "error_bound.h"
#include "fakeciphertext.h"
#include "eval_plain.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;

class DecompCache {
public:
    nlohmann::json cache_root;

    DecompCache(const std::string& cache_file) {
        cache_root = load_decomp_cache(cache_file);
        if (cache_root.is_null()) {
            throw std::runtime_error("Failed to load decomp cache: " + cache_file);
        }
    }

    EvalStep restore_EvalStep(const std::vector<double>& coeff) const {
        const std::string key = get_poly_type_key(coeff);
        if (!cache_root.contains(key)) {
            throw std::runtime_error("Cache key not found: " + key);
        }

        Poly poly(coeff);
        std::shared_ptr<Decomp> dcmp =
            reconstruct_decomp_from_cache(cache_root.at(key), poly, 0);

        if (!dcmp) {
            throw std::runtime_error("Failed to reconstruct Decomp from cache.");
        }

        return EvalStep(dcmp);
    }
};

EvalStep restore_EvalStep(const std::vector<double>& coeff, const std::string& cache_file) {
    DecompCache cache(cache_file);
    return cache.restore_EvalStep(coeff);
}

PYBIND11_MODULE(EBModule, m) {
    m.doc() = "pybind11 bindings for error_bound, parse, and fake ciphertext";

    py::class_<DecompCache>(m, "DecompCache")
        .def(py::init<const std::string&>(), py::arg("cache_file"))
        .def("restore_EvalStep", &DecompCache::restore_EvalStep, py::arg("coeff"));

    m.def("restore_EvalStep", &restore_EvalStep,
          py::arg("coeff"), py::arg("cache_file"));

    py::class_<fct::Ciphertext>(m, "Ciphertext")
        .def(py::init<>())
        .def(py::init<double, int, double, double, int>(),
             py::arg("x"), py::arg("size"), py::arg("err") = 0.0,
             py::arg("scale") = 1.0, py::arg("depth") = 0)
        .def(py::init<double, double, double, double, double, bool, int>(),
             py::arg("start"), py::arg("end"),
             py::arg("interval_size") = 1.0 / 16384.0, py::arg("err") = 0.0,
             py::arg("scale") = 1.0, py::arg("encode") = true, py::arg("depth") = 0)
        .def(py::init<const std::vector<std::pair<double, double>>&, double, double, bool, int>(),
             py::arg("intervals"), py::arg("err") = 0.0,
             py::arg("scale") = 1.0, py::arg("encode") = true, py::arg("depth") = 0)
        .def(py::init<const std::vector<double>&, double, double, bool, int>(),
             py::arg("x"), py::arg("err") = 0.0,
             py::arg("scale") = 1.0, py::arg("encode") = true, py::arg("depth") = 0)
        .def_readonly("g", &fct::Ciphertext::g)
        .def_readonly("r", &fct::Ciphertext::r)
        .def_readonly("scale", &fct::Ciphertext::scale)
        .def_readonly("depth", &fct::Ciphertext::depth)
        .def_property_readonly("intervals", [](const fct::Ciphertext& ct) {
            return ct.intervals ? *ct.intervals : std::vector<std::pair<double, double>>{};
        })
        .def("set_error", &fct::Ciphertext::set_error, py::arg("err"));

    m.def("evaluate", &fct::evaluate,
          py::arg("x"), py::arg("es"), py::arg("eb"));
    m.def("evaluate_composite", &fct::evaluate_composite,
          py::arg("x"), py::arg("eval_steps"), py::arg("eb"));
    m.def("estimate_range", &fct::estimate_range,
          py::arg("eb"), py::arg("interval"), py::arg("es"), py::arg("interval_size"));
    m.def("estimate_composite_range", &fct::estimate_composite_range,
          py::arg("eb"), py::arg("interval"), py::arg("eval_steps"), py::arg("interval_size"));
    m.def("generate_points", &fct::generate_points, py::arg("a"), py::arg("b"), py::arg("interval_size"));
    m.def("kappa", &fct::kappa, py::arg("coeff"), py::arg("interval"));
    m.def("max_abs", &fct::max_abs, py::arg("coeff"), py::arg("interval"));
    m.def("add", &fct::add, py::arg("ct1"), py::arg("ct2"), py::arg("eb"));
    m.def("add_pt", &fct::add_pt, py::arg("ct"), py::arg("scalar"));
    m.def("mult", &fct::mult, py::arg("ct1"), py::arg("ct2"), py::arg("eb"));
    m.def("square", &fct::square, py::arg("ct"), py::arg("eb"));
    m.def("mult_pt", &fct::mult_pt,
          py::arg("scalar"), py::arg("ct"), py::arg("eb"), py::arg("rescale") = true);
    m.def("dec", &fct::dec, py::arg("ct"));
    m.def("dec_interval", &fct::dec_interval, py::arg("ct"));
    m.def("dec_intervals", &fct::dec_intervals, py::arg("ct"));
          
    m.def("bit_diff", &bit_diff, py::arg("a"), py::arg("b"));
    m.def("clearNum", &clearNum, py::arg("a"));

    py::class_<ErrBound>(m, "ErrBound")
        .def(py::init<double, int, int, int>(),
             py::arg("sigma"), py::arg("N"), py::arg("h"), py::arg("s"))
        .def_readwrite("Bc", &ErrBound::Bc)
        .def_readwrite("Bs", &ErrBound::Bs)
        .def_readwrite("scale", &ErrBound::scale)
        .def("B_clean", &ErrBound::B_clean,
             py::arg("sigma"), py::arg("N"), py::arg("h"))
        .def("B_scale", &ErrBound::B_scale,
             py::arg("N"), py::arg("h"));

    py::class_<step>(m, "step")
        .def(py::init<>())
        .def_readwrite("op", &step::op)
        .def_readwrite("key1", &step::key1)
        .def_readwrite("key2", &step::key2)
        .def_readwrite("save_key", &step::save_key);

    py::class_<EvalStep>(m, "EvalStep")
        .def(py::init<>())
        .def_readwrite("term_count", &EvalStep::term_count)
        .def_readwrite("coeff_count", &EvalStep::coeff_count)
        .def_readwrite("powers", &EvalStep::powers)
        .def_readwrite("coeffs", &EvalStep::coeffs)
        .def_readwrite("terms", &EvalStep::terms)
        .def_readwrite("eval_step", &EvalStep::eval_step)
        .def("get_value", &EvalStep::get_value, py::arg("key"))
        .def("add_value", &EvalStep::add_value,
             py::arg("key_header"), py::arg("key_number"), py::arg("data"))
        .def("add_step", &EvalStep::add_step,
             py::arg("key1"), py::arg("key2"), py::arg("op"))
        .def("print_step", &EvalStep::print_step);
}
