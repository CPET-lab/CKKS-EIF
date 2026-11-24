"""Evaluate coefficient-type patterns by degree and optionally save OPD plans."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
from pathlib import Path
import sys

from tqdm import tqdm

REPO_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_DIR / "python"))

from OPD import NULL_DECOMP, cal_polyEval, make_all_poly_types
from OPD.src.save import get_poly_type_key, load_cache, save_cache, serialize_decomp

DEFAULT_WORKERS = 8
DEFAULT_BATCH_SIZE = 8
OPTIONS = ("all", "special")

_WORKER_INITIAL_CACHE = None
_WORKER_INTERMEDIATE_CACHE = None


def _init_worker() -> None:
    global _WORKER_INITIAL_CACHE, _WORKER_INTERMEDIATE_CACHE
    _WORKER_INITIAL_CACHE = {}
    _WORKER_INTERMEDIATE_CACHE = {}


def _evaluate_poly_type(coeff_type: list[str], divAssumption: bool) -> tuple[str, dict]:
    if _WORKER_INITIAL_CACHE is None or _WORKER_INTERMEDIATE_CACHE is None:
        _init_worker()
    result = cal_polyEval(
        coeff_type, divAssumption=divAssumption,
        initial_cache=_WORKER_INITIAL_CACHE, intermediate_cache=_WORKER_INTERMEDIATE_CACHE,
        cache_root_results=False,
    )
    key = get_poly_type_key(coeff_type)
    if result is NULL_DECOMP:
        raise RuntimeError(f"OPD failed for coefficient type {key}.")
    return key, serialize_decomp(result)


def _save_cache_atomic(cache: dict, output: Path) -> None:
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    try:
        save_cache(cache, str(temporary))
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()


def make_poly_types(degree: int, option: str) -> list[list[str]]:
    if degree < 1:
        raise ValueError("degree must be at least 1.")
    if option == "all":
        return make_all_poly_types(degree)
    if option != "special":
        raise ValueError(f"option must be one of {OPTIONS}.")
    return [
        ["F"] * (degree + 1),
        ["F" if exponent % 2 == degree % 2 else "0" for exponent in range(degree + 1)],
    ]


def default_output(option: str, divAssumption: bool=True) -> Path:
    if option not in OPTIONS:
        raise ValueError(f"option must be one of {OPTIONS}.")
    name = "decomp_cache" if option == "all" else "decomp_cache_special"
    if not divAssumption:
        name += "_no_div_assumption"
    return REPO_DIR / "cache" / "opd" / f"{name}.json"


def summarize_degree(degree: int, poly_types: list[list[str]], cache: dict) -> dict:
    distribution = Counter()
    examples = {}
    metrics = ("depth", "cmult", "pmult", "add")
    for coeff_type in poly_types:
        key = get_poly_type_key(coeff_type)
        complexity = cache[key]["complexity"]
        cost = tuple(complexity[name] for name in metrics)
        distribution[cost] += 1
        examples.setdefault(cost, key)
    return {
        "degree": degree,
        "total": len(poly_types),
        "ranges": {
            name: {"min": min(cost[i] for cost in distribution), "max": max(cost[i] for cost in distribution)}
            for i, name in enumerate(metrics)
        },
        "distribution": [
            dict(zip(metrics, cost), count=distribution[cost], example=examples[cost])
            for cost in sorted(distribution)
        ],
    }


def print_summary(summary: dict) -> None:
    print(f"\nDegree {summary['degree']}, total {summary['total']}")
    for name, values in summary["ranges"].items():
        print(f"{name}: {values['min']} ~ {values['max']}")
    print("\nDepth | CMult | PMult | Add | Count | Example")
    for row in summary["distribution"]:
        print(
            f"{row['depth']:>5} | {row['cmult']:>5} | {row['pmult']:>5} | "
            f"{row['add']:>3} | {row['count']:>5} | {row['example']}"
        )


def run_degrees(
    output: Path | None, min_degree: int=2, max_degree: int=8,
    divAssumption: bool=True, resume: bool=True, workers: int=DEFAULT_WORKERS,
    batch_size: int=DEFAULT_BATCH_SIZE, option: str="all",
) -> dict:
    if min_degree < 1 or max_degree < min_degree:
        raise ValueError("Use 1 <= min_degree <= max_degree.")
    if workers < 1 or batch_size < 1:
        raise ValueError("workers and batch_size must be at least 1.")
    if option not in OPTIONS:
        raise ValueError(f"option must be one of {OPTIONS}.")
    output = Path(output) if output is not None else None
    cache = load_cache(str(output)) if output is not None and resume else {}
    print(f"Option {option}: using {workers} worker processes with batch size {batch_size}.")

    for degree in range(min_degree, max_degree + 1):
        poly_types = make_poly_types(degree, option)
        pending = [coeff_type for coeff_type in poly_types if get_poly_type_key(coeff_type) not in cache]
        completed = len(poly_types) - len(pending)
        print(f"Degree {degree}: evaluating {len(pending)} new coefficient types ({completed} cached)...")
        with tqdm(total=len(poly_types), initial=completed, desc=f"Degree {degree}",
                  unit="case", mininterval=1.0, file=sys.stdout, ascii=True) as progress:
            for start in range(0, len(pending), batch_size):
                batch = pending[start:start + batch_size]
                before = completed
                try:
                    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker) as executor:
                        futures = {
                            executor.submit(_evaluate_poly_type, coeff_type, divAssumption): get_poly_type_key(coeff_type)
                            for coeff_type in batch
                        }
                        for future in as_completed(futures):
                            key, plan = future.result()
                            if key != futures[future]:
                                raise RuntimeError(f"Worker returned unexpected coefficient type {key}.")
                            cache[key] = plan
                            completed += 1
                            progress.update(1)
                except BaseException:
                    if output is not None and completed > before:
                        _save_cache_atomic(cache, output)
                    raise
                if output is not None:
                    _save_cache_atomic(cache, output)
        print_summary(summarize_degree(degree, poly_types, cache))
        if output is not None:
            print(f"Saved {completed}/{len(poly_types)} degree-{degree} cases to {output}")
    return cache


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, help="Evaluate one degree instead of a range.")
    parser.add_argument("--min-degree", type=int, help="First degree (default: 2).")
    parser.add_argument("--max-degree", type=int, help="Last degree (default: 8).")
    parser.add_argument("--save", action=argparse.BooleanOptionalAction, default=True,
                        help="Save plans as JSON; --no-save computes summaries without reading or writing a cache.")
    parser.add_argument("--option", choices=OPTIONS, default="all",
                        help="all: every coefficient type; special: all-F and matching-parity patterns.")
    parser.add_argument("--div-assumption", action=argparse.BooleanOptionalAction, default=True,
                        help="Treat nonzero non-leading quotient types as F (default: enabled).")
    parser.add_argument("--resume", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument("--output", type=Path, help="Output JSON path (default: an option-specific file under cache/opd/).")
    args = parser.parse_args()
    if args.degree is not None and (args.min_degree is not None or args.max_degree is not None):
        parser.error("Use --degree or --min-degree/--max-degree, not both.")
    min_degree = args.degree if args.degree is not None else (args.min_degree if args.min_degree is not None else 2)
    max_degree = args.degree if args.degree is not None else (args.max_degree if args.max_degree is not None else 8)
    output = (args.output if args.output is not None else default_output(args.option, args.div_assumption)) if args.save else None
    if min_degree < 1 or max_degree < min_degree:
        parser.error("Use 1 <= min_degree <= max_degree.")
    if args.workers < 1 or args.batch_size < 1:
        parser.error("workers and batch_size must be at least 1.")
    cache = run_degrees(
        output, min_degree=min_degree, max_degree=max_degree,
        divAssumption=args.div_assumption, resume=args.resume,
        workers=args.workers, batch_size=args.batch_size, option=args.option,
    )
    print(f"{'Saved' if args.save else 'Evaluated'} {len(cache)} total coefficient types.")


if __name__ == "__main__":
    main()
