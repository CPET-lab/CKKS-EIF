#!/usr/bin/env bash
set -euo pipefail

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository_root="$(cd "$script_directory/../.." && pwd)"
python_path="${PYTHON:-python3}"
range_figure_path="$script_directory/graph_ererange.pdf"
domain_figure_path="$script_directory/graph_eredomain.pdf"
# Used for automatic result selection: p = 2^p_exponent.
p_exponent=3
# Empty selects the latest completed combined-domain result for this p.
result_path=""
mode=both
layout=auto
# Plot only available FCT data for h=2^-start .. 2^-end (inclusive).
interval_size_start_exp=10
interval_size_end_exp=14

usage() {
    echo "Usage: $0 [RESULT_JSON [both|bounds|domain [auto|standard|test]]] [options]"
    echo "Generate graph_ererange.pdf and graph_eredomain.pdf in test/ERE/."
    echo "Defaults: the latest completed combined-domain result for p=2^$p_exponent and both final figures."
    echo "An explicit RESULT_JSON overrides automatic selection by p_exponent."
    echo "auto, standard, and test all use the final layout."
    echo "  --fct-first N  First plotted spacing exponent (default: $interval_size_start_exp)"
    echo "  --fct-last N   Last plotted spacing exponent (default: $interval_size_end_exp)"
}
position=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --help|-h) usage; exit 0 ;;
        --fct-first|--fct-last)
            if [[ $# -lt 2 ]]; then echo "$1 requires a value" >&2; exit 2; fi
            if [[ $1 == --fct-first ]]; then interval_size_start_exp="$2";
            else interval_size_end_exp="$2"; fi
            shift 2
            ;;
        --*) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
        *)
            case "$position" in
                0) result_path="$1" ;;
                1) mode="$1" ;;
                2) layout="$1" ;;
                *) usage >&2; exit 2 ;;
            esac
            position=$((position + 1))
            shift
            ;;
    esac
done
if [[ ! $interval_size_start_exp =~ ^[1-9][0-9]*$ || ! $interval_size_end_exp =~ ^[1-9][0-9]*$ ]]; then
    echo "FCT spacing exponents must be positive integers" >&2; exit 2
fi
if ((interval_size_start_exp > interval_size_end_exp)); then
    echo "FCT first must not exceed last" >&2; exit 2
fi
if [[ ! $p_exponent =~ ^[1-9][0-9]*$ ]] || ((p_exponent >= 31)); then
    echo "p_exponent must be a positive integer less than 31" >&2; exit 2
fi
export MPLCONFIGDIR="${MPLCONFIGDIR:-$script_directory/result/.matplotlib}"
mkdir -p "$MPLCONFIGDIR"
if [[ -z $result_path ]]; then
    result_path="$("$python_path" -B - "$script_directory/src" "$script_directory/result" "$p_exponent" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
from plot_range import latest_result
print(latest_result(sys.argv[2], p_exponent=int(sys.argv[3])))
PY
)"
fi
if [[ $result_path != /* ]]; then result_path="$PWD/$result_path"; fi
if [[ ! -f $result_path ]]; then echo "Result JSON not found: $result_path" >&2; exit 2; fi
case "$mode" in both|bounds|domain) ;; *) usage >&2; exit 2 ;; esac
case "$layout" in auto|standard|test) ;; *) usage >&2; exit 2 ;; esac
bounds_plotter="$script_directory/src/plot_range.py"
echo "Input: $result_path"
if [[ $mode == both || $mode == bounds ]]; then
    "$python_path" -B "$bounds_plotter" "$result_path" \
        "$range_figure_path" \
        --fct-first "$interval_size_start_exp" --fct-last "$interval_size_end_exp"
fi
if [[ $mode == both || $mode == domain ]]; then
    "$python_path" -B "$script_directory/src/plot_domain.py" "$result_path" \
        "$domain_figure_path"
fi
