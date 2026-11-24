# OPD implementation

`src/algorithm.py` implements OPD. The remaining modules define polynomial and plan
types, power generation, decomposition formatting, and JSON cache handling.

With the repository's `python/` directory on the module search path, import the
public API from `OPD`, for example `from OPD import Poly, cal_polyEval`.

Executable experiments and their usage guide are in
[`test/OPD/readme.md`](../../test/OPD/readme.md).
