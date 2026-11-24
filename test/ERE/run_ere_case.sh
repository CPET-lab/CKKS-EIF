#!/usr/bin/env bash
set -euo pipefail

# Edit these defaults, then run this script without arguments.
script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository_root="$(cd "$script_directory/../.." && pwd)"
p_exponent=3
rounds=10
workers=1
domain_mode=combined
ct_spacing_exp=16
interval_size_start_exp=10
interval_size_end_exp=14
coefficient_directory="$repository_root/cache/ops"
result_directory="$script_directory/result"
executable_path="$repository_root/build/EREExperiment"
ct_input_path=""
output_name=""
quick=false
check_only=false

usage() {
    echo "Usage: $0 [P_EXPONENT [ROUNDS [WORKERS]]] [options]"
    echo "Generate CT/FCT/Baseline JSON in result/; no plots are generated."
    echo "Default: combined domain, fresh CT, FCT h=2^-10..2^-14."
    echo "  --ct-input PATH   Reuse saved CT; recompute Baseline/FCT without CT measurement"
    echo "  --fresh-ct        Measure new CT on the selected evaluation domain"
    echo "  --domain-mode M   combined (default) or split input regions"
    echo "  --fct-first N     First FCT spacing exponent (default: $interval_size_start_exp)"
    echo "  --fct-last N      Last FCT spacing exponent (default: $interval_size_end_exp)"
    echo "  --ct-exp N        CT spacing exponent (default: $ct_spacing_exp)"
    echo "  --quick           One CT round at spacing 2^-12; write *_quick.json"
    echo "  --output-name NAME Save JSON under this filename in result/"
    echo "  --check-only      Validate paths/settings without building or running"
}

position=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --help|-h) usage; exit 0 ;;
        --ct-input|--fct-first|--fct-last|--ct-exp|--output-name|--domain-mode)
            if [[ $# -lt 2 ]]; then echo "$1 requires a value" >&2; exit 2; fi
            case "$1" in
                --ct-input) ct_input_path="$2" ;;
                --fct-first) interval_size_start_exp="$2" ;;
                --fct-last) interval_size_end_exp="$2" ;;
                --ct-exp) ct_spacing_exp="$2" ;;
                --output-name) output_name="$2" ;;
                --domain-mode) domain_mode="$2" ;;
            esac
            shift 2
            ;;
        --quick) quick=true; shift ;;
        --fresh-ct) ct_input_path=""; shift ;;
        --check-only) check_only=true; shift ;;
        --*) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
        *)
            case "$position" in
                0) p_exponent="$1" ;;
                1) rounds="$1" ;;
                2) workers="$1" ;;
                *) usage >&2; exit 2 ;;
            esac
            position=$((position + 1))
            shift
            ;;
    esac
done

for value in "$p_exponent" "$rounds" "$workers" "$ct_spacing_exp" \
    "$interval_size_start_exp" "$interval_size_end_exp"; do
    if [[ ! $value =~ ^[1-9][0-9]*$ ]]; then
        echo "p, rounds, workers, and spacing exponents must be positive integers" >&2
        exit 2
    fi
done
if [[ $domain_mode != combined && $domain_mode != split ]]; then
    echo "Domain mode must be combined or split" >&2; exit 2
fi
if ((p_exponent >= 31 || interval_size_start_exp > interval_size_end_exp)); then
    echo "Require p exponent < 31 and FCT first <= last" >&2; exit 2
fi
coefficient_path="$coefficient_directory/coeff_${p_exponent}_30_depth.txt"
for path in "$coefficient_path" "$repository_root/cache/opd/decomp_cache.json"; do
    if [[ ! -f $path ]]; then echo "Input file not found: $path" >&2; exit 2; fi
done
if [[ $quick == true && $ct_input_path == auto ]]; then ct_input_path=""; fi
if [[ $ct_input_path == auto ]]; then
    ct_suffix=""
    if [[ $p_exponent -ne 2 ]]; then ct_suffix="_p${p_exponent}"; fi
    ct_input_path="$result_directory/ere_results${ct_suffix}.json"
fi
if [[ -n $ct_input_path ]]; then
    if [[ $ct_input_path != /* ]]; then ct_input_path="$PWD/$ct_input_path"; fi
    if [[ ! -f $ct_input_path ]]; then echo "CT input JSON not found: $ct_input_path" >&2; exit 2; fi
    if [[ $quick == true ]]; then echo "--quick cannot be combined with --ct-input" >&2; exit 2; fi
fi
if [[ $quick == true ]]; then rounds=1; ct_spacing_exp=12; fi
suffix=""
if [[ $p_exponent -ne 2 ]]; then suffix="_p${p_exponent}"; fi
if [[ $domain_mode == combined ]]; then suffix="_combined${suffix}";
elif [[ -n $ct_input_path ]]; then suffix="${suffix}_test"; fi
if [[ $quick == true ]]; then suffix="${suffix}_quick"; fi
output_name="${output_name:-ere_results${suffix}.json}"
if [[ ! $output_name =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*\.json$ ]]; then
    echo "--output-name must be a JSON filename without directories" >&2; exit 2
fi
result_path="$result_directory/$output_name"
if [[ -n $ct_input_path && $(realpath -m "$ct_input_path") == $(realpath -m "$result_path") ]]; then
    echo "CT input and output must be different files" >&2; exit 2
fi
if [[ -n $ct_input_path ]]; then
    echo "ERE p=2^$p_exponent domain=$domain_mode CT=reused($ct_input_path) FCT=2^-$interval_size_start_exp..2^-$interval_size_end_exp"
else
    echo "ERE p=2^$p_exponent domain=$domain_mode rounds=$rounds workers=$workers CT=2^-$ct_spacing_exp FCT=2^-$interval_size_start_exp..2^-$interval_size_end_exp"
fi
echo "Output: $result_path"
if [[ $check_only == true ]]; then exit 0; fi

cd "$repository_root"
cmake -S "$repository_root" -B "$repository_root/build" -DCMAKE_BUILD_TYPE=Release
cmake --build "$repository_root/build" --config Release --target EREExperiment --parallel 2
mkdir -p "$result_directory"
arguments=(--domain-mode "$domain_mode" --p-exp "$p_exponent" --rounds "$rounds"
    --workers "$workers" --ct-exp "$ct_spacing_exp"
    --fct-first "$interval_size_start_exp" --fct-last "$interval_size_end_exp"
    --coeff-dir "$coefficient_directory" --output "$result_path")
if [[ -n $ct_input_path ]]; then arguments+=(--ct-input "$ct_input_path"); fi
"$executable_path" "${arguments[@]}"
