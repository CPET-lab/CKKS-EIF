# ERE guide

## What ERE computes

Error-aware Range Estimation (ERE) estimates the output interval of a polynomial
evaluated over CKKS, accounting for the error model adopted in the paper.
Its inputs are an evaluation sequence, an input interval, CKKS error parameters,
and a subinterval width `h`.

ERE follows the evaluation sequence using Fake-Ciphertexts, each carrying a
nominal polynomial and an error radius. It partitions the input interval,
bounds polynomial variation on each subinterval, and expands the bounds by the
propagated error radius. The output is an interval enclosing the modeled CKKS
output; reducing `h` reduces the interval-approximation overestimation.

The experiment compares three ranges and saves them as JSON:

| Result | Meaning |
| --- | --- |
| CT | Observed range from sampled actual CKKS evaluations |
| FCT | ERE-estimated range for each selected `h` |
| Baseline | Plaintext extrema expanded once by the rescaling-error bound |

## Running the test programs

Use the environment described in the [root README](../../README.md).
Run commands from the repository root (`/app` inside Docker).

Before running an ERE experiment, `cache/opd/decomp_cache.json` must contain
an OPD evaluation plan matching the coefficient-type pattern of every polynomial
stage to be evaluated. Generate any missing plans with `save_OPD.py` as described
in the [OPD guide](../OPD/readme.md). The supplied cache covers all degree 2–8
coefficient-type patterns.

### run_ere_case.sh: run the experiment

Check input paths and settings without running the experiment:

```bash
bash test/ERE/run_ere_case.sh --check-only
```

Run with the default settings:

```bash
bash test/ERE/run_ere_case.sh
```

Defaults are `p=8`, 10 CT rounds, one worker, CT spacing `2^-16`, and FCT widths
`2^-10` through `2^-14`. For a shorter run, use `--quick` (one CT round with
CT spacing `2^-12`).

To change the p exponent, CT rounds, workers, and FCT widths:

```bash
bash test/ERE/run_ere_case.sh 2 10 1 --fct-first 10 --fct-last 14
```

The positional arguments are `P_EXPONENT ROUNDS WORKERS`, with `p=2^P_EXPONENT`.
Inputs are `cache/opd/decomp_cache.json` and
`cache/ops/coeff_<P_EXPONENT>_30_depth.txt`.
The default result is `test/ERE/result/ere_results_combined_p3.json`.
Console `lower=PASS` and `upper=PASS` indicate that the FCT bounds enclose
the observed CT range on the corresponding sides.

Use `--ct-input PATH` to reuse saved CT results while recomputing FCT and
Baseline. The launcher supports `--help`.

### plot_ere_case.sh: plot the results

Plot the latest completed combined-domain result for the configured p
(default: `p=8`):

```bash
PYTHON=python3 bash test/ERE/plot_ere_case.sh
```

To select a result file and plot type explicitly:

```bash
PYTHON=python3 bash test/ERE/plot_ere_case.sh \
  test/ERE/result/ere_results_combined_p3.json bounds
```

Choose `bounds`, `domain`, or `both` (default). The script saves
`graph_ererange.pdf` and `graph_eredomain.pdf`, with PNG copies, under
`test/ERE/` according to the selected plot type. Use `--help` for all options.
