"""Run OPS using the paper settings and release input caches."""

import argparse
import os
from pathlib import Path
import sys
from types import SimpleNamespace

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "python"))

from common.paths import OPS_RESULTS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--criteria", choices=("depth", "time"), default="depth")
    parser.add_argument("--pruning-type", choices=("base", "prune", "lex"), default="base")
    parser.add_argument("--start", type=int, default=8)
    parser.add_argument("--end", type=int, default=8)
    parser.add_argument("--e-num", type=int, default=30)
    parser.add_argument("--eb-params", nargs=4,
                        metavar=("SIGMA", "N", "HAMMING_WEIGHT", "SCALE_BITS"),
                        default=(3.2, 2 ** 17, 128, 50))
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--max-n", type=int, default=8)
    parser.add_argument("--max-search-depth", type=int, default=25)
    parser.add_argument("--interval-exp", type=int, default=14)
    parser.add_argument("--save-log", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--save-result", action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args()
    try:
        args.eb_params = (float(args.eb_params[0]), int(args.eb_params[1]),
                          int(args.eb_params[2]), int(args.eb_params[3]))
    except ValueError:
        parser.error("--eb-params requires SIGMA followed by integer N, HAMMING_WEIGHT, SCALE_BITS")
    if args.start < 1 or args.start > args.end:
        parser.error("Use 1 <= start <= end")
    if min(args.workers, args.max_search_depth, args.interval_exp) < 1 or args.max_n < 2:
        parser.error("workers, max-search-depth, and interval-exp must be positive; max-n >= 2")
    available = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    options = SimpleNamespace(
        **vars(args), prune_range=("c", "p", "a"), multiThread=True,
        cpu_count=min(args.workers, available), pruning_min_routes=10_000,
        interval_size=2 ** -args.interval_exp,
    )
    OPS_RESULTS.mkdir(parents=True, exist_ok=True)
    if options.save_log:
        suffix = f"_{options.pruning_type}" if options.criteria == "depth" else ""
        options.FILENAME = str(OPS_RESULTS / f"output_{options.e_num}_{options.criteria}{suffix}.txt")
    else:
        options.FILENAME = None
    from OPS.src.search_sep import multi_remez
    multi_remez(options=options)


if __name__ == "__main__":
    main()
