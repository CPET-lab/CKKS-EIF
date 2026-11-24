"""Shared visual and y-axis styles for the ERE paper figures."""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
from matplotlib.ticker import Formatter, MaxNLocator
import numpy as np


DOMAIN_TITLES = {
    "target": r"$\mathcal{I}_{\mathrm{test}}$",
    "non_target": r"$\mathcal{Z}_{\mathrm{test}}$",
}
TITLE_SIZE = 10
LABEL_SIZE = 10
TICK_SIZE = 9
ANNOTATION_SIZE = 8
GRID_ALPHA = 0.25
CURVE_COLOR = "tab:blue"


def _exponent(values: np.ndarray) -> int:
    magnitude = float(np.max(np.abs(values)))
    return math.floor(math.log10(magnitude)) if magnitude else 0


def _auto_significant_digits(values: np.ndarray) -> int:
    magnitude = float(np.max(np.abs(values)))
    span = float(np.ptp(values))
    if not magnitude or not span:
        return 3
    return max(3, min(10, math.ceil(math.log10(magnitude / span)) + 2))


def _scientific_label(value: float, significant_digits: int) -> str:
    coefficient, exponent = f"{abs(value):.{significant_digits - 1}e}".split("e")
    coefficient = coefficient.rstrip("0").rstrip(".")
    return f"{coefficient}e{int(exponent)}"


class _ScaledFormatter(Formatter):
    def __init__(self, exponent: int, significant_digits: int) -> None:
        self.exponent = exponent
        self.significant_digits = significant_digits

    def __call__(self, value: float, pos: int | None = None) -> str:
        return f"{value / 10**self.exponent:.{self.significant_digits}g}"

    def get_offset(self) -> str:
        return f"1e{self.exponent}" if self.exponent else ""


class _NearOneFormatter(Formatter):
    def __init__(self, significant_digits: int, values_are_offsets: bool) -> None:
        self.significant_digits = significant_digits
        self.values_are_offsets = values_are_offsets

    def __call__(self, value: float, pos: int | None = None) -> str:
        delta = value if self.values_are_offsets else value - 1.0
        if delta == 0:
            return "1"
        sign = "+" if delta > 0 else "-"
        return f"1{sign}{_scientific_label(delta, self.significant_digits)}"


class _ScientificFormatter(Formatter):
    def __init__(self, significant_digits: int) -> None:
        self.significant_digits = significant_digits

    def __call__(self, value: float, pos: int | None = None) -> str:
        if value == 0:
            return "0"
        sign = "-" if value < 0 else ""
        return f"{sign}{_scientific_label(value, self.significant_digits)}"


def style_near_one_y_axis(
    axis: plt.Axes,
    *,
    significant_digits: int = 3,
    values_are_offsets: bool = False,
    max_ticks: int = 4,
) -> None:
    """Show ticks near one as 1±ceN; offsets may be plotted directly."""
    if significant_digits < 1:
        raise ValueError("significant_digits must be positive")
    axis.yaxis.set_major_locator(MaxNLocator(nbins=max_ticks))
    axis.yaxis.set_major_formatter(
        _NearOneFormatter(significant_digits, values_are_offsets)
    )
    axis.tick_params(axis="y", labelsize=TICK_SIZE)


def style_scaled_y_axis(
    axis: plt.Axes,
    reference_values: np.ndarray,
    *,
    significant_digits: int | None = None,
    max_ticks: int = 4,
) -> None:
    """Show ticks as coefficients with one shared 1eN factor."""
    values = np.asarray(reference_values, dtype=float)
    if values.size == 0:
        raise ValueError("reference_values must not be empty")
    if significant_digits is None:
        significant_digits = _auto_significant_digits(values)
    if significant_digits < 1:
        raise ValueError("significant_digits must be positive")
    axis.yaxis.set_major_locator(MaxNLocator(nbins=max_ticks))
    axis.yaxis.set_major_formatter(
        _ScaledFormatter(_exponent(values), significant_digits)
    )
    axis.tick_params(axis="y", labelsize=TICK_SIZE)


def style_scientific_y_axis(
    axis: plt.Axes,
    reference_values: np.ndarray,
    *,
    significant_digits: int | None = None,
    max_ticks: int = 4,
) -> None:
    """Write the scientific exponent on each tick, without an axis factor."""
    values = np.asarray(reference_values, dtype=float)
    if values.size == 0:
        raise ValueError("reference_values must not be empty")
    if significant_digits is None:
        significant_digits = _auto_significant_digits(values)
    if significant_digits < 1:
        raise ValueError("significant_digits must be positive")
    axis.yaxis.set_major_locator(MaxNLocator(nbins=max_ticks))
    axis.yaxis.set_major_formatter(_ScientificFormatter(significant_digits))
    axis.tick_params(axis="y", labelsize=TICK_SIZE)


def style_domain_axis(axis: plt.Axes, domain: str, *, show_title: bool = True) -> None:
    if show_title:
        axis.set_title(DOMAIN_TITLES[domain], fontsize=TITLE_SIZE)
    axis.tick_params(axis="x", labelsize=TICK_SIZE)
    axis.grid(True, axis="y", alpha=GRID_ALPHA)
