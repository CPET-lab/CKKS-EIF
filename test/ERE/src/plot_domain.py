"""Final Figure 1: the composed function on the combined input domain."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext, MaxNLocator
import mpmath as mp
import numpy as np
from plot_style import CURVE_COLOR, GRID_ALPHA, LABEL_SIZE, TICK_SIZE

EXPERIMENT_DIR = Path(__file__).resolve().parent.parent
DEFAULT_RESULT = EXPERIMENT_DIR / "result/ere_results_combined.json"
DEFAULT_FIGURE = EXPERIMENT_DIR / "graph_eredomain.pdf"


def load_result(result_path):
    data = json.loads(Path(result_path).read_text(), parse_float=str)
    if not data.get("case", {}).get("stages") or not data.get("parameters", {}).get("evaluation_domains"):
        raise ValueError("Polynomial stages and evaluation domains are required")
    return data


def evaluate_composition(x, stages):
    values = np.asarray(x, dtype=float)
    for stage in stages:
        values = np.polynomial.polynomial.polyval(values, np.asarray(stage["coefficients"], dtype=float))
    return values


def target_value_offset(x: np.ndarray, stages: list[dict]) -> np.ndarray:
    # Evaluate F(x)-1 before converting to float to preserve its variation.
    with mp.workdps(80):
        coefficients = [
            [mp.mpf(value) for value in stage["coefficients"]]
            for stage in stages
        ]

        def composed(value: mp.mpf) -> mp.mpf:
            for stage_coefficients in coefficients:
                value = mp.polyval(stage_coefficients[::-1], value)
            return value

        return np.asarray([
            float(composed(mp.mpf(str(value))) - 1)
            for value in x
        ])


def plot(data, output, linthresh=1e-10):
    output = Path(output)
    end = float(data["parameters"]["p"]) - 1 + float(data["parameters"]["clean_error"])
    x = np.linspace(0, end, 131073)
    y = evaluate_composition(x, data["case"]["stages"])
    figure, axis = plt.subplots(figsize=(4.25, 2.7))
    axis.set_yscale("symlog", base=10, linthresh=linthresh, linscale=0.7)
    axis.plot(x, y, color=CURVE_COLOR, linewidth=2)
    axis.axhline(0, color="0.25", linewidth=1, linestyle="--", zorder=1)
    axis.set(xlim=(0, end), ylim=(-4e-5, 2), xlabel="$x$", ylabel="$F(x)$")
    if data["parameters"]["p"] == 4:
        axis.set_xticks(np.arange(0, 3.01, 0.5))
    else:
        axis.xaxis.set_major_locator(MaxNLocator(nbins=7, integer=True))
    axis.set_yticks([-1e-5, -1e-8, 0, 1e-8, 1e-5, 1e-2, 1])
    axis.yaxis.set_major_formatter(LogFormatterMathtext(base=10, linthresh=linthresh))
    axis.tick_params(labelsize=TICK_SIZE)
    axis.xaxis.label.set_fontsize(LABEL_SIZE)
    axis.yaxis.label.set_fontsize(LABEL_SIZE)
    axis.grid(which="major", alpha=GRID_ALPHA)
    axis.spines[["top", "right"]].set_visible(False)
    figure.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.98)
    output.parent.mkdir(parents=True, exist_ok=True)
    for extension in [".pdf", ".png"]:
        figure.savefig(output.with_suffix(extension), bbox_inches="tight", dpi=200)
    plt.close(figure)


def plot_composite_function(data, output=DEFAULT_FIGURE):
    plot(data, output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", nargs="?", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("output", nargs="?", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--output", dest="output_override", type=Path)
    parser.add_argument("--linthresh", type=float, default=1e-10)
    args = parser.parse_args()
    if not np.isfinite(args.linthresh) or args.linthresh <= 0:
        parser.error("linthresh must be finite and positive")
    plot(load_result(args.source or args.result), args.output_override or args.output, args.linthresh)
    print(f"Saved {(args.output_override or args.output).with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
