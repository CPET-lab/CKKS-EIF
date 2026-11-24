# CKKS-EIF

Research code for **encrypted indicator function (EIF) evaluation over CKKS**,
based on the paper:

> **Efficient Indicator Function Evaluation over CKKS via Cost- and Error-Aware Polynomial Composition**
>
> Junhyuk Kim and Dongwoo Kim, Dongguk University.

Paper: [TBD](TBD)

## Protocols

- **OPD — Optimal Polynomial Decomposition:** Finds low-cost polynomial
  evaluation plans through decomposition.
  See the [OPD guide](test/OPD/readme.md).
- **ERE — Error-aware Range Estimation:** Estimates polynomial output ranges
  while accounting for CKKS evaluation errors.
  See the [ERE guide](test/ERE/readme.md).
- **OPS — Optimal Path Search:** Searches for polynomial compositions that meet
  a target precision using cost and range estimates from OPD and ERE.
  See the [OPS guide](test/OPS/readme.md).

## Repository layout

| Path | Contents |
| --- | --- |
| `python/OPD` | OPD implementation |
| `python/OPS` | OPS implementation |
| `python/common` | Shared paths, runtime-cache loading, and Remez approximation |
| `cpp/ERE` | ERE implementation and Python bindings |
| `cpp/common` | Shared polynomial and CKKS implementation |
| `cache/opd` | Polynomial decomposition plans |
| `cache/runtime` | Runtime estimates |
| `cache/ops` | Selected polynomial coefficients |
| `test/OPD`, `test/ERE`, `test/OPS` | Experiment entry points |

## Environment

Docker provides the Python environment, modified SEAL library, and C++ build
tools. Install [Docker Engine or Docker Desktop](https://docs.docker.com/get-started/get-docker/)
on the host and confirm that `docker version` can contact the Docker server.

| Component | Configuration |
| --- | --- |
| Verified platform | Linux/amd64 |
| Container OS | Ubuntu 24.04 |
| Python | 3.12 |
| C++ | C++20 |
| CMake | 3.18 or newer |
| SEAL | CPET_SEAL_Light, based on SEAL 4.1.2 |
| Python packages | Versions listed in [requirements.txt](requirements.txt) |

The Dockerfile pins [CPET_SEAL_Light](https://github.com/CPET-lab/CPET_SEAL_Light)
to commit `a4a428c0d630a31a2e4179fd964019f2e3b3264d`, with Gaussian noise and
Intel HEXL disabled. OpenMP, OpenBLAS, and MKL threads are set to one per process.

The commands below use a POSIX shell, such as Bash on Linux, macOS, or WSL2.
Windows Docker Desktop must use Linux containers. Network access is required
to clone the repository and download dependencies during image construction.

## Installation

Clone the repository and build the Docker image:

```bash
git clone https://github.com/CPET-lab/CKKS-EIF.git
cd CKKS-EIF
docker build --progress=plain -t ckks-eif:local .
```

The default build uses two compiler processes. To use one compiler process:

```bash
docker build --build-arg BUILD_JOBS=1 --progress=plain -t ckks-eif:local .
```

`BUILD_JOBS` controls compilation; experiment parallelism is configured
separately through each program's `--workers` option.

Open a shell in the image:

```bash
docker run --rm -it ckks-eif:local
```

The container's working directory is `/app`. Python, SEAL, and the compiler are
installed inside the image. Experiment commands are described in the protocol
guides linked above.

## Acknowledgments

The multi-interval Remez algorithm in `python/common/remez` was independently
implemented based on the method described by
[Lee et al. (EUROCRYPT 2021)](https://doi.org/10.1007/978-3-030-77870-5_22).

This work was supported by Institute of Information & Communications Technology
Planning & Evaluation (IITP) under the Artificial Intelligence Convergence
Innovation Human Resources Development (IITP-2025-RS-2023-00254592) grant funded by
the Korean government (MSIT).

## Citation

If you use this code in your research, please cite the accompanying paper:

```bibtex
@misc{kim_ckks_eif,
  author = {Kim, Junhyuk and Kim, Dongwoo},
  title = {Efficient Indicator Function Evaluation over CKKS via Cost- and Error-Aware Polynomial Composition}
}
```

<!-- Add the publication year, venue, paper URL, and DOI once finalized. -->

## License

Project code is licensed under the [MIT License](LICENSE).
The bundled nlohmann/json header retains its MIT license and copyright notices;
see [THIRD_PARTY_NOTICES](THIRD_PARTY_NOTICES). Additional dependencies retain
their own licenses.

## Contact

Junhyuk Kim: [jhkim@dgu.ac.kr](mailto:jhkim@dgu.ac.kr).
