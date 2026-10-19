"""Replot the four Kim/Hornung convergence figures at manuscript text size.

Example from the thesis-workspace root::

    python3 thesis/04_OpenFOAM_HiSA_VnV/Figures/plot_shock_convergence.py \
        --data-root '../../Files/ValData/FIN' \
        --audit-dir tmp/shock_caption_font_review

The 170 mm fixed PDF width equals the manuscript linewidth. Text is 9.5--10
points at that width, approximately the caption size. No bounding-box scaling
is applied. Original pressure exports, source scripts and PNGs are untouched.
The source readers/detectors are imported directly, including the distinct
Hornung 8M HPC detector and the original history-marker selections. Missing
data cause the affected figure to be skipped, never reconstructed from images.
Requires numpy, pandas, matplotlib and the Latin Modern OpenType fonts.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, ticker
import numpy as np
import pandas as pd


WIDTH_MM = 170.0
BLUE = "#1f77b4"


def source_module(directory, filename, module_name):
    """Load helpers without running any source script's plotting entry point."""
    spec = importlib.util.spec_from_file_location(module_name, directory / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def configure_fonts(directory):
    for name in ("regular", "italic", "bold", "bolditalic"):
        font_manager.fontManager.addfont(directory / f"lmroman10-{name}.otf")
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["Latin Modern Roman"],
        "mathtext.fontset": "cm", "text.usetex": False,
        "font.size": 9.5, "axes.labelsize": 10, "axes.titlesize": 10,
        "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
        "legend.fontsize": 9.5, "axes.linewidth": 0.7,
        "pdf.fonttype": 3,  # Embed native CFF OTF outlines, as in the profile panels.
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })


def save_figure(fig, output, audit_dir):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    # A font increase must not clip labels at the fixed media-box boundary.
    for ax in fig.axes:
        box = ax.get_tightbbox(renderer)
        if (box.x0 < -0.5 or box.y0 < -0.5 or
                box.x1 > fig.bbox.width + 0.5 or box.y1 > fig.bbox.height + 0.5):
            raise ValueError(f"Labels exceed the fixed figure boundary: {output.name}")
    fig.savefig(output, metadata={"Creator": Path(__file__).name,
                                 "Title": output.stem})
    if audit_dir:
        fig.savefig(audit_dir / (output.stem + ".png"), dpi=200)
    plt.close(fig)


def load_kim_mesh(directory, source):
    rows = []
    for name, cells in zip(["6k", "16k", "38k", "90k", "250k", "750k", "2p3M"],
                           [6000, 16000, 38000, 90000, 250000, 750000, 2300000]):
        profiles = sorted((directory / name).glob("*k.csv"), key=source.extract_iteration)
        if not profiles:
            raise FileNotFoundError(f"No original profile in {directory / name}")
        path = profiles[-1]
        x, p = source.read_stagline_csv(path)
        result = source.find_shock_position((x - source.X_NOSE) / source.D, p / source.P_INF)
        if result is None:
            raise ValueError(f"No shock detected in {path}")
        position, spacing = result
        rows.append({"case": name, "cells": cells, "source": str(path),
                     "distance_mm": -position * source.D * 1000,
                     "half_interval_mm": spacing * source.D * 500,
                     "relative_error_percent": (-position - 0.228) / 0.228 * 100,
                     "half_interval_percent": spacing / 2 / 0.228 * 100})
    return pd.DataFrame(rows)


def plot_mesh(data, kind):
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH_MM / 25.4, 64 / 25.4))
    fig.subplots_adjust(left=0.09, right=0.987, bottom=0.205, top=0.97, wspace=0.35)
    kim = kind == "Kim"
    cells = data.cells.to_numpy() if kim else data.equivalent_cells.to_numpy()
    distance = data.distance_mm if kim else data.delta_mm
    half = data.half_interval_mm if kim else data.half_sampling_interval_mm
    relative_half = data.half_interval_percent if kim else data.half_sampling_interval_percent
    reference = 2.280 if kim else 12.189
    specs = [
        (distance, half, reference, "^", r"$d$ [mm]" if kim else r"$\Delta$ [mm]",
         (1.95, 2.40) if kim else (9, 12.5), "%.2f" if kim else "%.1f",
         f"{kind} ref. ({reference:.3f} mm)"),
        (data.relative_error_percent, relative_half, 0, "o", "Relative error [%]",
         (-5, 5) if kim else (-25, 2.5), "%.0f", f"{kind} ref. (0 %)"),
    ]
    for ax, (y, error, ref, marker, label, limits, fmt, legend) in zip(axes, specs):
        ax.set_xscale("log")
        ax.set_xlim((4500, 3e6) if kim else (11000, 1.1e7))
        ax.set_xticks(cells)
        ax.set_xticklabels([f"{c/1e6:g}M" if c >= 1e6 else f"{c/1e3:g}k" for c in cells])
        ax.xaxis.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
        ax.xaxis.set_minor_formatter(ticker.NullFormatter())
        ax.axhline(ref, color=BLUE, ls="--", lw=1.1, label=legend, zorder=2)
        points = ax.errorbar(cells, y, yerr=error, fmt=f"-{marker}", color=BLUE,
                            ms=4.5, capsize=2.5, lw=1.1, elinewidth=0.8,
                            label="CFD (present)", zorder=3)
        assert np.array_equal(points.lines[0].get_ydata(), np.asarray(y))
        ax.set_ylim(*limits)
        # Keep the original data ranges, with fewer labelled ticks for readability.
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5 if kim else 7))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter(fmt))
        ax.set_ylabel(label, labelpad=3)
        ax.set_xlabel("Number of cells" if kim else "Equivalent number of cells", labelpad=3)
        ax.tick_params(which="major", direction="in", length=3.5, width=0.7, top=True, right=True)
        ax.tick_params(which="minor", direction="in", length=2, width=0.5, top=True, right=True)
        ax.grid(True, which="major", color="0.85", lw=0.5)
        ax.set_axisbelow(True)
        leg = ax.legend(loc="lower right", framealpha=1, fancybox=False, edgecolor="0.6",
                        borderpad=0.35, labelspacing=0.3, handlelength=1.8, handletextpad=0.5)
        leg.get_frame().set_linewidth(0.5)
    return fig


def history_axes(ax, title, ylabel, hornung):
    ax.set_title(title, pad=5)
    ax.set_xlabel(r"Pseudo-iterations [$10^3$]" if hornung else r"Iteration [$10^3$]", labelpad=3)
    ax.set_ylabel(ylabel, labelpad=3)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
    formatter = ticker.ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)
    ax.yaxis.set_major_formatter(formatter)
    ax.grid(True, color="0.86", linewidth=0.5)
    ax.set_axisbelow(True)
    if hornung:
        ax.tick_params(direction="in", top=True, right=True, length=3, width=0.7)
    ax.margins(x=0.07, y=0.17)


def plot_kim_history(data):
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH_MM / 25.4, 64 / 25.4))
    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.205, top=0.89, wspace=0.35)
    for ax, mesh, label, color in zip(axes, ["250k", "2p3M"],
                                    ["(a) Kim, 250k cells", "(b) Kim, 2.3M cells"],
                                    ["#ff0000", "#9467bd"]):
        history = data[data.mesh == mesh]
        ax.errorbar(history.iteration / 1000, history.d_mm,
                    yerr=history.half_sampling_interval_mm, fmt="o-", color=color,
                    ms=3.5, capsize=2, linewidth=1.1, elinewidth=0.8, zorder=3)
        ax.axhline(history.iloc[-1].d_mm, color="0.45", ls="--", lw=0.8, zorder=1)
        history_axes(ax, label, r"$d$ [mm]", False)
    return fig


def plot_hornung_history(data, source):
    fig, axes = plt.subplots(2, 2, figsize=(WIDTH_MM / 25.4, 108 / 25.4))
    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.12, top=0.93,
                        wspace=0.35, hspace=0.57)
    for index, (ax, (name, label, cells)) in enumerate(zip(axes.flat, source.PLOT_CASES)):
        history = data[data.case == name]
        x, y = history.iteration_thousands, history.shock_standoff_mm
        show = source.marker_mask(history, name)
        ax.axhline(y.iloc[-1], color="0.4", ls="--", lw=0.8, zorder=2)
        line, = ax.plot(x, y, color=BLUE, lw=1.1, zorder=3)
        ax.errorbar(x[show], y[show], yerr=history.half_sampling_interval_mm[show],
                    fmt="o", color=BLUE, ms=3.5, elinewidth=0.8, capsize=2, zorder=4)
        assert np.array_equal(line.get_ydata(), y.to_numpy())
        history_axes(ax, f"({chr(97+index)}) {label}", r"$\Delta$ [mm]", True)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--font-dir", type=Path,
                        default=Path("/mnt/c/texlive/2026/texmf-dist/fonts/opentype/public/lm"))
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--audit-dir", type=Path)
    args = parser.parse_args()
    configure_fonts(args.font_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.audit_dir:
        args.audit_dir.mkdir(parents=True, exist_ok=True)
    report = {"width_mm": WIDTH_MM, "text_size_pt": [9.5, 10], "figures": {}}
    for kind, directory_name in [("Kim", "Kim/01 mesh_sens"), ("Hornung", "Hornung/mesh_sens")]:
        directory = args.data_root.resolve() / directory_name
        for plot_type in ["mesh", "iteration"]:
            name = f"03_{kind}_{plot_type}_convergence"
            try:
                # The original scripts use the same helper-module name; replace it
                # for each family to avoid accidentally mixing geometry/detectors.
                helper = source_module(directory, "plot_stagline_pressure_iters.py",
                                       "plot_stagline_pressure_iters")
                if plot_type == "mesh":
                    if kind == "Kim":
                        data = load_kim_mesh(directory, helper)
                    else:
                        source = source_module(directory, "plot_hornung_convergence.py", "hornung_mesh")
                        data = source.load_results()
                    fig = plot_mesh(data, kind)
                else:
                    source = source_module(directory, f"plot_{kind.lower()}_iteration_convergence.py",
                                           f"{kind.lower()}_history")
                    data = source.load_histories()
                    fig = plot_kim_history(data) if kind == "Kim" else plot_hornung_history(data, source)
                save_figure(fig, args.output_dir / (name + ".pdf"), args.audit_dir)
                if args.audit_dir:
                    data.to_csv(args.audit_dir / (name + "_data.csv"), index=False)
                report["figures"][name] = {"status": "replotted", "source_directory": str(directory)}
                print(f"Replotted {name}: original data, 170 mm width, 9.5--10 pt text")
            except (FileNotFoundError, ValueError) as error:
                plt.close("all")
                report["figures"][name] = {"status": "skipped", "reason": str(error)}
                print(f"SKIPPED {name}: {error}", file=sys.stderr)
    if args.audit_dir:
        (args.audit_dir / "replot_report.json").write_text(json.dumps(report, indent=2) + "\n")
    if any(item["status"] == "skipped" for item in report["figures"].values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
