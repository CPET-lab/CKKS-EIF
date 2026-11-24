"""Print an OPD decomposition for numerical coefficients or coefficient types."""

import argparse
import ast
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "python"))

from OPD import NULL_DECOMP, Poly, cal_polyEval
from OPD.src.print import type_decomp_to_str


def parse_coefficients(raw: str) -> list[float]:
    text = raw.strip()
    if text.startswith("["):
        try:
            values = ast.literal_eval(text)
        except (SyntaxError, ValueError) as error:
            raise ValueError("Expected a list of finite numerical coefficients.") from error
    else:
        values = [float(value.strip()) for value in text.split(",")]
    if not isinstance(values, list) or not values:
        raise ValueError("At least one coefficient is required.")
    if any(type(value) not in (int, float) for value in values):
        raise ValueError("Coefficients must be numbers, not strings or booleans.")
    coefficients = [float(value) for value in values]
    if not all(math.isfinite(value) for value in coefficients):
        raise ValueError("Coefficients must be finite.")
    if coefficients[-1] == 0:
        raise ValueError("The last coefficient must be nonzero; trailing zeros are not accepted.")
    return coefficients


def parse_coefficient_types(raw: str) -> list[str]:
    text = raw.strip()
    if text.startswith("["):
        try:
            values = ast.literal_eval(text)
        except (SyntaxError, ValueError) as error:
            raise ValueError("Expected a list containing only '0', 'I', and 'F'.") from error
    else:
        values = text.split(",") if "," in text else list(text)
    if not isinstance(values, list) or not values:
        raise ValueError("At least one coefficient type is required.")
    if any(not isinstance(value, str) for value in values):
        raise ValueError("Coefficient types must be strings: '0', 'I', or 'F'.")
    types = [value.strip().upper() for value in values]
    if any(value not in {"0", "I", "F"} for value in types):
        raise ValueError("Coefficient types must contain only 0, I, and F.")
    if types[-1] == "0":
        raise ValueError("The last coefficient type must be I or F; trailing zeros are not accepted.")
    return types


def calculate_optimal(poly_input):
    result = cal_polyEval(
        Poly(poly_input), divAssumption=True,
        cache_root_results=False,
    )
    if result is NULL_DECOMP:
        raise RuntimeError("No decomposition was found.")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--coeff", help="Coefficients from constant to leading term: '[1, 2.5, 0, 4]' or 1,2.5,0,4")
    inputs.add_argument("--type", help="Coefficient types: '[\"I\", \"F\", \"0\", \"I\"]', I,F,0,I, or IF0I")
    args = parser.parse_args()
    try:
        values = parse_coefficients(args.coeff) if args.coeff is not None else parse_coefficient_types(args.type)
        result = calculate_optimal(values)
    except (ValueError, TypeError, OverflowError, RuntimeError) as error:
        parser.error(str(error))
    print(f"Input: {values}")
    if args.type is not None:
        symbols = ", ".join(f"a{i}: {kind}" for i, kind in enumerate(values) if kind != "0")
        print(f"Coefficient types (I = nonzero integer, F = non-integer): {symbols or 'zero polynomial'}")
        print(f"OPD decomposition: {type_decomp_to_str(result)}")
    else:
        print(f"OPD decomposition: {result.restore_dcmp()}")
    print(f"Complexity (Depth|CMult|PMult|Add): {result.comp.return_params()}")


if __name__ == "__main__":
    main()
