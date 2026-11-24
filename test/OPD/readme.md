# OPD guide

## What OPD computes

Optimal Polynomial Decomposition (OPD) searches for a low-cost evaluation plan
for a polynomial under the CKKS cost model. The input can be numerical
coefficients or a coefficient-type pattern: `0` (zero), `I` (nonzero integer),
and `F` (non-integer), ordered from the constant term to the leading term.

OPD explores recursive decompositions of the form `f(x) = tau(x)*p(x) + q(x)`,
where the factor is `x^i` or a coefficient-scaled power `a*x^i`. It tracks
previously computed powers and their depths so that later operations can reuse
them, and selects the minimum-cost plan among the candidates it explores.

The output is a polynomial decomposition and its cost tuple,
`Depth|CMult|PMult|Add`:

| Component | Meaning |
| --- | --- |
| Depth | Multiplicative depth |
| CMult | Ciphertext-ciphertext multiplications |
| PMult | Non-integer plaintext scalar multiplications |
| Add | Additions |

Costs are compared lexicographically in this order: Depth, CMult, PMult, Add.
Type inputs produce symbolic decompositions; numerical inputs produce
decompositions using the supplied coefficients. Batch execution also reports
cost ranges and distributions, and can save evaluation plans as JSON.

For type input, `divAssumption=True` treats nonzero non-leading quotient
coefficients in an `a*x^i` decomposition as F. Numerical input always uses
the actual floating-point division results.

## Running the test programs

Use Python 3.12 and run commands from the repository root (`/app` inside Docker).
`save_OPD.py` also requires `tqdm`, included in `requirements.txt`.
See the [root README](../../README.md) for environment setup.

### test_OPD.py: one polynomial

Choose numerical coefficients (`--coeff`) or coefficient types (`--type`):

```bash
python3 test/OPD/test_OPD.py --coeff '[1, 2.5, 0, 4]'
python3 test/OPD/test_OPD.py --type IF0I
```

Both inputs have the pattern of `1 + 2.5*x + 4*x^3` and return cost `2|2|0|2`.
The program prints the decomposition and its cost without reading or updating
a saved cache. Inputs must be nonempty and end in a nonzero coefficient or type;
numerical coefficients must be finite.

### save_OPD.py: evaluate and save by degree

Evaluate all coefficient-type patterns of one degree without reading or writing
a saved cache:

```bash
python3 test/OPD/save_OPD.py --degree 4 --no-save
```

Evaluate and save all patterns for degrees 2 through 8:

```bash
python3 test/OPD/save_OPD.py --min-degree 2 --max-degree 8
```

Each degree `n` has `2 * 3^n` patterns (19,674 total for degrees 2–8).
The program prints cost ranges and distributions. Plans are saved to
`cache/opd/decomp_cache.json` after each batch; existing entries are skipped
by default. The supplied cache already contains the degree 2–8 results.

Use `--output PATH` for another file, `--no-resume` to recompute, or
`--workers N` to set parallelism (default: 8). `--no-div-assumption` considers
both I and F for quotient coefficients and uses a separate default cache.

Each program supports `--help`.
