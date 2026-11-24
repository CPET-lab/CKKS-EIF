cmake_minimum_required(VERSION 3.18)
project(ApproxRelease LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Threads REQUIRED)
find_package(SEAL 4.1 CONFIG REQUIRED)
find_package(Python 3.12 COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

add_library(ere_core STATIC
    cpp/ERE/error_bound.cpp
    cpp/ERE/fakeciphertext.cpp
    cpp/common/eval_plain.cpp
    cpp/common/parse.cpp
    cpp/common/polyEval_class.cpp)
set_target_properties(ere_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(ere_core PUBLIC cpp/ERE cpp/common)

add_library(ckks_core STATIC
    cpp/common/CKKS_params.cpp
    cpp/common/evaluate.cpp
    cpp/common/multicore.cpp)
target_link_libraries(ckks_core PUBLIC ere_core SEAL::seal Threads::Threads)

pybind11_add_module(EBModule NO_EXTRAS cpp/ERE/bindings.cpp)
target_link_libraries(EBModule PRIVATE ere_core)

add_executable(EREExperiment test/ERE/src/ere_experiment.cpp)
target_link_libraries(EREExperiment PRIVATE ckks_core)
target_compile_definitions(EREExperiment PRIVATE APPROX_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
