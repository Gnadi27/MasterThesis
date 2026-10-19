"""Rebuild the two full-width profile panels of the Kim extraction figure.

The original profiles and shock-location calculations are retained. Only the
figure geometry, legend placement and font sizes change. The authoritative
CSV reader and shock detector are loaded directly from the original script at
``DATA_DIR/plot_stagline_pressure_iters.py``; no duplicate implementation is
maintained here.

Example from the thesis-workspace root (Windows-hosted data through WSL)::

    python3 thesis/04_OpenFOAM_HiSA_VnV/Figures/plot_kim_shock_extraction.py \
        --data-dir '../../Files/ValData/FIN/Kim/01 mesh_sens' \
        --font-dir /mnt/c/texlive/2026/texmf-dist/fonts/opentype/public/lm

Requires numpy, pandas and matplotlib. Uses Latin Modern OpenType fonts and
Computer Modern mathtext directly, without a running LaTeX installation.
The original PNG assets and simulation data are never modified.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import runpy

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
import numpy as np


WIDTH_MM = 170.0
HEIGHT_MM = 46.0
CASES = (
    ("6k", "25k.csv", "coarse"),
    ("2p3M", "45k.csv", "fine"),
)


def configure_fonts(font_dir: Path) -> None:
    """Embed the thesis's actual Latin Modern text fonts in each PDF."""
    for name in ("lmroman10-regular.otf", "lmroman10-italic.otf",
                 "lmroman10-bold.otf", "lmroman10-bolditalic.otf"):
        font_manager.fontManager.addfont(font_dir / name)
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Latin Modern Roman"],
        "font.size": 9.5,
        "axes.labelsize": 10,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "text.usetex": False,
        "mathtext.fontset": "cm",
        "pdf.fonttype": 3,  # Native CFF OTF glyph outlines; avoids invalid Type 42 wrapping.
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
    })


def add_compact_inset(ax, source, x_nd, p_nd, x_shock, dx_local, color):
    """Retain the source inset samples, limits and markings in a wider box."""
    x_mm = (source["X_NOSE"] + source["D"] * x_nd) * 1000.0
    x_shock_mm = (source["X_NOSE"] + source["D"] * x_shock) * 1000.0
    dx_mm = dx_local * source["D"] * 1000.0
    half_span_mm = 8.0 * dx_mm
    local = np.abs(x_mm - x_shock_mm) <= half_span_mm

    inset = ax.inset_axes([0.66, 0.24, 0.30, 0.51])
    inset.plot(x_mm[local], p_nd[local], color=color, alpha=0.85,
               linewidth=1.1, marker="o", markersize=3.0)
    inset.axvspan(x_shock_mm - dx_mm / 2.0, x_shock_mm + dx_mm / 2.0,
                  color="red", alpha=0.15, zorder=0)
    inset.axvline(x_shock_mm, color="red", linestyle=":", linewidth=1.4)
    inset.set_xlim(x_shock_mm - half_span_mm, x_shock_mm + half_span_mm)
    inset.set_ylim(0.0, np.ceil(np.max(p_nd[local]) / 5.0) * 5.0)
    inset.xaxis.set_major_locator(MaxNLocator(nbins=3))
    inset.xaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    inset.set_yticks([0, 20, 40])
    inset.set_xlabel(r"$x$ [mm]", fontsize=8, labelpad=1,
                     bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.3})
    inset.set_ylabel(r"$p/p_\infty$", fontsize=8, labelpad=1)
    inset.set_title("Shock detail", fontsize=8, pad=2)
    inset.tick_params(labelsize=8, length=2, pad=2)
    inset.grid(True, alpha=0.3)


def plot_panel(data_dir: Path, output_dir: Path, source, case, csv_name, suffix):
    csv_path = data_dir / case / csv_name
    # These are the unchanged authoritative source functions.
    x, p = source["read_stagline_csv"](csv_path)
    x_nd = (x - source["X_NOSE"]) / source["D"]
    p_nd = p / source["P_INF"]
    result = source["find_shock_position"](x_nd, p_nd)
    if result is None:
        raise ValueError(f"No shock could be detected in {csv_path}")
    x_shock, dx_local = result
    x_shock_mm = (source["X_NOSE"] + source["D"] * x_shock) * 1000.0
    half_mm = dx_local * source["D"] * 1000.0 / 2.0
    mask = (x_nd >= source["X_MIN_D"]) & (x_nd <= source["X_MAX_D"])

    # Each source case has one iteration file: Normalize(vmin=vmax) maps it
    # to turbo(0), preserving the original dark-purple profile color exactly.
    color = plt.get_cmap("turbo")(0.0)
    fig, ax = plt.subplots(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    fig.subplots_adjust(left=0.078, right=0.985, bottom=0.265, top=0.815)
    (line,) = ax.plot(x[mask] * 1000.0, p_nd[mask], linewidth=1.3,
                     alpha=0.85, color=color, label="present CFD result")
    ax.axvspan(x_shock_mm - half_mm, x_shock_mm + half_mm,
               color="red", alpha=0.15, zorder=0)
    ax.axvline(x_shock_mm, color="red", linestyle=":", linewidth=1.4,
               label=rf"detected normal shock $({x_shock_mm:.3f} \pm {half_mm:.3f})\,\mathrm{{mm}}$")
    if case == "2p3M":
        add_compact_inset(ax, source, x_nd, p_nd, x_shock, dx_local, color)
    ax.axvline(source["X_NOSE"] * 1000.0, color="black", linestyle="--",
               linewidth=1.0, label=r"$x_\mathrm{nose}$")
    ax.set_xlabel(r"$x$ [mm]", labelpad=2)
    ax.set_ylabel(r"$p/p_\infty$", labelpad=2)
    ax.set_yticks([0, 10, 20, 30, 40, 50])
    ax.tick_params(length=3, pad=2)
    ax.grid(True, alpha=0.3)
    fig.legend(*ax.get_legend_handles_labels(), loc="upper center",
               bbox_to_anchor=(0.54, 1.0), ncol=3, fontsize=9.0,
               frameon=False, handlelength=2.0, handletextpad=0.5,
               columnspacing=1.1, borderaxespad=0.2)

    # Assert that plotting did not resample, smooth, or alter the source data.
    assert np.array_equal(line.get_xdata(), x[mask] * 1000.0)
    assert np.array_equal(line.get_ydata(), p_nd[mask])
    out_path = output_dir / f"03_Kim_stag_line_shock_identified_{suffix}_wide.pdf"
    # Fixed media box: no tight bounding-box crop or later aspect stretching.
    fig.savefig(out_path, metadata={"Title": f"Kim cylinder {suffix} mesh pressure profile",
                                   "Creator": Path(__file__).name})
    plt.close(fig)
    print(f"{out_path.name}: {WIDTH_MM:g} x {HEIGHT_MM:g} mm; "
          f"{int(mask.sum())} plotted samples; x_shock={x_shock_mm:.9f} mm; "
          f"half interval={half_mm:.9f} mm; d/D={-x_shock:.9f}; "
          "profile coordinate differences from source=0")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True,
                        help="Original Kim mesh_sens directory, including its plotting script.")
    parser.add_argument("--font-dir", type=Path,
                        default=Path("/mnt/c/texlive/2026/texmf-dist/fonts/opentype/public/lm"),
                        help="Directory containing the Latin Modern Roman OpenType fonts.")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent,
                        help="Destination for the two new wide PDF assets.")
    args = parser.parse_args()
    data_dir = args.data_dir.resolve()
    source = runpy.run_path(str(data_dir / "plot_stagline_pressure_iters.py"),
                           run_name="kim_pressure_profile_helpers")
    configure_fonts(args.font_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for case, csv_name, suffix in CASES:
        plot_panel(data_dir, args.output_dir, source, case, csv_name, suffix)


if __name__ == "__main__":
    main()
