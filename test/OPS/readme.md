# OPS guide

## What OPS computes

Optimal Path Search (OPS) searches for a composition of polynomials that
approximates an indicator function to the target precision under the paper's
CKKS error model. It extends a composition one polynomial at a time, using
multi-interval Remez approximation to map the current target and non-target
intervals toward 1 and 0.

OPD supplies each polynomial's evaluation plan and operation cost. ERE estimates
the resulting output intervals, which become the input intervals of the next
stage. OPS selects a composition satisfying the target precision according to
the requested optimization criterion:

| `--criteria` | Objective |
| --- | --- |
| `depth` | Lexicographic `Depth\|CMult\|PMult\|Add` cost |
| `time` | CKKS evaluation time estimated from the runtime cache |

The result is an ordered sequence of polynomials, reported through their degrees,
saved coefficients, total operation cost, and estimated CKKS evaluation time.
The search log also records the measured OPS search time.

## Running the test programs

Use the environment and built `EBModule` described in the
[root README](../../README.md). Run commands from the repository root
(`/app` inside Docker).

Before searching, `cache/opd/decomp_cache.json` must contain plans matching the
candidate coefficient-type patterns. The supplied cache covers all degree 2–8
types; see the [OPD guide](../OPD/readme.md) to generate missing plans.
`cache/runtime/time_cache.json` must also cover the requested depths, degrees,
and polynomial types. OPS uses the [ERE implementation](../ERE/readme.md)
to estimate output ranges.

### test.py: search for a polynomial composition

Run a search for `p=4` with either optimization criterion:

```bash
python3 test/OPS/test.py --criteria depth --start 2 --end 2 --workers 2
python3 test/OPS/test.py --criteria time --start 2 --end 2 --workers 2
```

`--start` and `--end` are inclusive exponents of p (`p=2^k`). Defaults are
exponents 8 through 8, target precision `2^-30`, maximum component degree 8,
maximum total depth 25, ERE width `2^-14`, and one worker.
`--max-n`, `--max-search-depth`, and `--interval-exp` change these bounds.

The default CKKS parameters are sigma=3.2, N=131072, secret-key Hamming weight
128, and scale `2^50`. For the 15-bit setting with scale `2^35`, use:

```bash
python3 test/OPS/test.py --criteria depth --start 2 --end 6 --e-num 15 \
  --eb-params 3.2 131072 128 35
```

Selected coefficients are saved to `cache/ops/coeff_<k>_<e>_<criterion>.txt`,
where `k` is the p exponent and `e` is the target precision exponent.
A completed run replaces an existing file with the same name.
Add `--no-save-result` to preserve the supplied coefficient files.

Logs are saved to `test/OPS/output_<e>_depth_<mode>.txt` or
`test/OPS/output_<e>_time.txt`. Depth modes are selected with
`--pruning-type base|prune|lex` (default: `base`). Use `--no-save-log` for
terminal output. Run one search at a time per release directory because
coefficient files, logs, and temporary route storage are shared.

Use `--help` to display all search options.
