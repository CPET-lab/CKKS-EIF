"""Figure 2: bounds with a separated upper scale and a fitted lower scale."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
from matplotlib.ticker import FuncFormatter, ScalarFormatter
import numpy as np
from plot_style import ANNOTATION_SIZE, CURVE_COLOR, GRID_ALPHA, LABEL_SIZE, TICK_SIZE, TITLE_SIZE

EXPERIMENT_DIR = Path(__file__).resolve().parent.parent
DEFAULT_RESULT = EXPERIMENT_DIR / "result/ere_results_combined.json"
DEFAULT_FIGURE = EXPERIMENT_DIR / "graph_ererange.pdf"


def _load_result(path):
    data = json.loads(Path(path).read_text())
    parameters = data.get("parameters", {})
    if "ere_option" in parameters:
        raise ValueError("Regenerate this result with the single ERE implementation")
    p = parameters.get("p")
    if not isinstance(p, int) or p < 2 or len(parameters.get("evaluation_domains", [])) != 1:
        raise ValueError("Figure 2 requires combined-domain results with p >= 2; run run_ere_case.sh")
    case = data.get("case", {})
    if not case.get("fct_ranges") or "baseline_range" not in case or "ct_range" not in case:
        raise ValueError("Baseline, CT, and FCT ranges are required")
    return data


def latest_result(directory, *, p_exponent=None):
    candidates = []
    for path in Path(directory).glob("*.json"):
        try:
            data = _load_result(path)
        except (OSError, ValueError):
            continue
        parameters = data["parameters"]
        if p_exponent is None or parameters["p"] == 2 ** p_exponent:
            candidates.append(path)
    if not candidates:
        selection = "" if p_exponent is None else f" for p=2^{p_exponent}"
        raise FileNotFoundError(f"No completed combined-domain result{selection} found in {directory}")
    return max(candidates, key=lambda path: (path.stat().st_mtime_ns, path.name))


def _select_fct_intervals(data, first_exponent=None, last_exponent=None):
    if any(value is not None and value < 1 for value in (first_exponent, last_exponent)):
        raise ValueError("FCT spacing exponents must be positive")
    if first_exponent is not None and last_exponent is not None and first_exponent > last_exponent:
        raise ValueError("FCT first must not exceed last")
    rows = sorted((row for row in data["case"]["fct_ranges"]
                   if (first_exponent is None or row["exponent"] >= first_exponent)
                   and (last_exponent is None or row["exponent"] <= last_exponent)),
                  key=lambda row: row["exponent"])
    if not rows:
        raise ValueError("No FCT data exists in the requested spacing range")
    return {**data, "case": {**data["case"], "fct_ranges": rows}}


def save(figure, output):
    figure.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(output.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.close(figure)


def axis_break(axis, start, end):
    axis.fill_between([0, 1], start, end, transform=axis.transAxes,
                      color="white", linewidth=0, zorder=5)
    phase = np.linspace(-np.pi, np.pi, 65)
    center = (start + end) / 2
    for side in [0, 1]:
        x = side + 0.006 * np.sin(2 * phase)
        y = center + np.linspace(-0.025, 0.025, len(phase))
        axis.plot(x, y, transform=axis.transAxes, color="white", linewidth=5,
                  clip_on=False, zorder=8)
        axis.plot(x, y, transform=axis.transAxes, color="black", linewidth=1.1,
                  clip_on=False, zorder=9)


def upper_tick(value, _):
    if value == 0:
        return "$1$"
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = abs(value) / 10.0 ** exponent
    coefficient = "" if np.isclose(mantissa, 1.0) else rf"{mantissa:.3g}\times"
    sign = "+" if value > 0 else "-"
    return rf"$1{sign}{coefficient}10^{{{exponent}}}$"


def upper_broken_axis(axis, reference, fct):
    if np.min(reference) <= 0 or np.min(fct) <= 4 * np.max(reference):
        return False
    lower_low = 10.0 ** np.floor(np.log10(np.min(reference)))
    lower_high = 10.0 ** np.ceil(np.log10(np.max(reference)))
    if lower_high <= lower_low:
        lower_high = lower_low * 10
    upper_low = np.min(fct) / 1.3
    if upper_low <= lower_high:
        return False
    magnitude = 10.0 ** np.floor(np.log10(np.max(fct)))
    upper_high = np.ceil(1.2 * np.max(fct) / magnitude) * magnitude
    lower_end, upper_start = 0.30, 0.40
    lower_log_span = np.log10(lower_high / lower_low)
    log_span = np.log10(upper_high / upper_low)

    # Both retained bands are logarithmic; only the empty gap is compressed.
    def forward(values):
        values = np.asarray(values)
        lower = lower_end * np.log10(
            np.maximum(values, np.finfo(float).tiny) / lower_low) / lower_log_span
        gap = lower_end + (upper_start - lower_end) * (values - lower_high) / (upper_low - lower_high)
        upper = upper_start + (1 - upper_start) * np.log10(
            np.maximum(values, np.finfo(float).tiny) / upper_low) / log_span
        return np.where(values <= lower_high, lower, np.where(values < upper_low, gap, upper))

    def inverse(values):
        values = np.asarray(values)
        lower = lower_low * 10.0 ** (values * lower_log_span / lower_end)
        gap = lower_high + (values - lower_end) * (upper_low - lower_high) / (upper_start - lower_end)
        upper = upper_low * 10.0 ** ((values - upper_start) * log_span / (1 - upper_start))
        return np.where(values <= lower_end, lower, np.where(values < upper_start, gap, upper))

    axis.set_yscale("function", functions=(forward, inverse))
    axis.set_ylim(lower_low, upper_high)
    reference_ticks = 10.0 ** np.arange(np.log10(lower_low), np.log10(lower_high) + 1)
    ticks = 10.0 ** np.arange(np.ceil(np.log10(upper_low)), np.floor(np.log10(upper_high)) + 1)
    axis.set_yticks(np.unique(np.r_[reference_ticks, ticks]))
    axis.yaxis.set_major_formatter(FuncFormatter(upper_tick))
    axis_break(axis, lower_end, upper_start)
    return True


def lower_broken_axis(axis, reference, fct):
    values = np.unique(np.r_[reference, fct])
    if len(values) < 5:
        return False
    gaps = np.diff(values)
    cuts = np.sort(np.argsort(gaps)[-2:])
    other_gaps = np.delete(gaps, cuts)
    reference_gap = np.ptp(reference)
    if reference_gap <= 0 or np.min(gaps[cuts]) <= 4 * max(np.max(other_gaps), reference_gap):
        return False
    groups = np.split(values, cuts + 1)
    if not (groups[-1][0] <= np.min(reference) <= np.max(reference) <= groups[-1][-1]):
        return False
    if reference_gap < 0.05 * np.ptp(groups[-1]):
        return False
    padding = max(0.3 * reference_gap, 0.04 * np.ptp(groups[-1]))
    bounds = np.array([[group[0] - padding, group[-1] + padding] for group in groups])
    if np.any(bounds[:-1, 1] >= bounds[1:, 0]):
        return False
    knots = bounds.ravel()
    positions = np.array([0.0, 0.12, 0.18, 0.30, 0.36, 1.0])

    def interpolate(values, source, target):
        values = np.asarray(values)
        indices = np.clip(np.searchsorted(source, values, side="right") - 1, 0, len(source) - 2)
        fraction = (values - source[indices]) / (source[indices + 1] - source[indices])
        return target[indices] + fraction * (target[indices + 1] - target[indices])

    axis.set_yscale("function", functions=(
        lambda values: interpolate(values, knots, positions),
        lambda values: interpolate(values, positions, knots),
    ))
    axis.set_ylim(knots[0], knots[-1])
    axis.set_yticks(np.r_[groups[0].mean(), groups[1].mean(),
                         np.linspace(groups[-1][0], groups[-1][-1], 3)])
    for start, end in [(0.12, 0.18), (0.30, 0.36)]:
        axis_break(axis, start, end)
    return True


def bound_zoom(axis, baseline, ct, styles, boundary, source_x, rectangle=None):
    if rectangle is None:
        rectangle = [0.055, 0.30, 0.53, 0.28] if boundary == 1 else [0.36, 0.43, 0.49, 0.25]
    inset = axis.inset_axes(rectangle, zorder=10)
    offset = 1.0 if boundary == 1 else 0.0
    values = np.array([baseline, ct]) - offset
    low, high = sorted(values)
    delta = abs(ct - baseline)
    padding = max(0.45 * delta, 8 * abs(np.spacing(max(abs(low), abs(high)))))
    inset.set_ylim(low - padding, high + padding)
    inset.set_xlim(0, 1)
    for (_, color, linestyle), value in zip(styles[:2], values):
        inset.axhline(value, color=color, linestyle=linestyle, linewidth=1.8)
    inset.set_xticks([])
    inset.set_yticks([])
    inset.annotate("", xy=(0.15, high), xytext=(0.15, low),
                   fontsize=ANNOTATION_SIZE,
                   arrowprops={"arrowstyle": "<->", "color": "0.25", "linewidth": 0.7,
                               "mutation_scale": 4, "shrinkA": 0, "shrinkB": 0})
    coefficient, exponent = f"{delta:.2e}".split("e")
    inset.text(0.28, low + delta / 2,
               rf"$\Delta={coefficient}\times10^{{{int(exponent)}}}$",
               va="center", ha="left", fontsize=ANNOTATION_SIZE)
    axis.add_artist(ConnectionPatch(
        xyA=(0.5, 0.0 if boundary == 1 else 1.0), coordsA=inset.transAxes,
        xyB=(source_x, low + delta / 2), coordsB=axis.transData,
        color="0.5", linewidth=1.0, linestyle="--", zorder=1,
    ))


def linear_bounds_axis(axis, values):
    low, high = float(np.min(values)), float(np.max(values))
    padding = max(0.08 * (high - low),
                  np.finfo(float).eps * max(abs(low), abs(high), 1e-12) * 8)
    axis.set_yscale("linear")
    axis.set_ylim(low - padding, high + padding)
    axis.set_yticks(np.linspace(low - padding, high + padding, 5))
    formatter = ScalarFormatter(useOffset=True, useMathText=True)
    formatter.set_powerlimits((0, 0))
    axis.yaxis.set_major_formatter(formatter)


def plot_bounds(data, output):
    output = Path(output)
    data = data.get("case", data)
    exponents = np.array([row["exponent"] for row in data["fct_ranges"]])
    fct = np.array([row["range"] for row in data["fct_ranges"]])
    baseline = np.asarray(data["baseline_range"])
    ct = np.asarray(data.get("ct_reference_range", data.get("ct_range")))
    if not np.all(np.isfinite(np.r_[baseline, ct, fct.ravel()])):
        raise ValueError("Bounds must be finite")
    lower_values = np.r_[baseline[0], ct[0], fct[:, 0]]
    figure, axes = plt.subplots(2, 1, figsize=(4.25, 4.2), sharex=True)
    styles = [
        ("Baseline", CURVE_COLOR, "--"),
        ("CT", "tab:orange", "-"),
        ("FCT", "tab:green", "-"),
    ]
    for axis, boundary, title in zip(axes, [1, 0], ["Upper bound", "Lower bound"]):
        offset = 1.0 if boundary == 1 else 0.0
        for (label, color, linestyle), values in zip(
            styles, [np.full(len(exponents), baseline[boundary]),
                     np.full(len(exponents), ct[boundary]), fct[:, boundary]],
        ):
            axis.plot(exponents, values - offset, label=label, color=color,
                      linestyle=linestyle, marker="^" if label == "FCT" else None,
                      linewidth=1.8, markersize=6,
                      zorder=6 if label == "FCT" else 2)
        axis.set_title(title, fontsize=TITLE_SIZE)
        axis.set_ylabel("Output bound", fontsize=LABEL_SIZE)
        broken_upper = False
        broken_lower = False
        if boundary == 1:
            reference = np.array([baseline[1], ct[1]]) - offset
            upper_values = np.r_[baseline[1], ct[1], fct[:, 1]] - offset
            broken_upper = upper_broken_axis(axis, reference, fct[:, 1] - offset)
            if not broken_upper:
                if np.min(upper_values) > 0:
                    axis.set_yscale("log")
                    first = np.floor(np.log10(np.min(upper_values)))
                    last = max(first + 1, np.ceil(np.log10(np.max(upper_values))))
                    axis.set_ylim(10.0 ** first, 10.0 ** last)
                    axis.set_yticks(10.0 ** np.arange(first, last + 1))
                else:
                    linear_bounds_axis(axis, upper_values)
                axis.yaxis.set_major_formatter(FuncFormatter(upper_tick))
        else:
            reference = np.array([baseline[0], ct[0]])
            linear_bounds_axis(axis, lower_values)
            magnitude = np.max(np.abs(lower_values))
            exponent = int(np.floor(np.log10(magnitude))) if magnitude else 0
            unit = 10.0 ** exponent
            axis.yaxis.set_major_formatter(FuncFormatter(
                lambda value, _, unit=unit: rf"${value / unit:.6f}$"
            ))
            axis.text(0, 1.02, rf"$\times10^{{{exponent}}}$", transform=axis.transAxes,
                      fontsize=ANNOTATION_SIZE)
        axis.minorticks_off()
        axis.grid(axis="y", alpha=GRID_ALPHA)
        axis.tick_params(labelsize=TICK_SIZE)
        baseline_fails = (baseline[1] < ct[1] if boundary == 1
                          else baseline[0] > ct[0])
        if baseline_fails:
            axis.annotate("Baseline FAIL", xy=(0.98, baseline[boundary] - offset),
                          xycoords=axis.get_yaxis_transform(),
                          xytext=(0, 12 if boundary == 1 else 4),
                          textcoords="offset points", ha="right", va="bottom",
                          color=CURVE_COLOR, fontsize=ANNOTATION_SIZE, zorder=10)
        bound_zoom(axis, baseline[boundary], ct[boundary], styles, boundary,
                   exponents[-1] - 0.6 if broken_lower else exponents[0] + 0.75,
                   rectangle=[0.055, 0.64, 0.53, 0.28] if broken_lower else None)
    axes[-1].set_xticks(exponents, [rf"$2^{{-{value}}}$" for value in exponents])
    axes[-1].set_xlabel("$h$", fontsize=LABEL_SIZE)
    axes[-1].set_xlim(exponents[0] - 0.15, exponents[-1] + 0.15)
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="upper center", ncol=3, fontsize=9)
    figure.tight_layout(rect=(0, 0, 1, 0.93), h_pad=1.4, pad=0.6)
    save(figure, output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", nargs="?", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("output", nargs="?", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--fct-first", type=int, default=8)
    parser.add_argument("--fct-last", type=int, default=12)
    args = parser.parse_args()
    data = _select_fct_intervals(_load_result(args.source or args.result), args.fct_first, args.fct_last)
    plot_bounds(data, args.output)
    print(f"Saved {args.output.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
