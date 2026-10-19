#!/usr/bin/env python3
"""Render Figure 7.22 with 10-point text at the thesis's 170 mm print width.

The adjacent CSV retains the original 3% cubic finite-stencil table, unchanged.
Original WSL source paths and hashes are recorded in provenance.json. Run this
script with Python, matplotlib, numpy and pandas; outputs are written one level
above this directory. The PDF uses fixed dimensions (no tight-bounding-box
rescaling), so width=\\linewidth preserves the intended printed font size.
"""
from pathlib import Path
import argparse
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, ticker
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
COLORS = {"T": "#0072B2", "U": "#D07800",
          "muEff": "#7650B4", "kappaEff": "#BE558E"}
LABELS = {"T": r"$T$ gradient (DNS)", "U": r"$U$ gradient (DNS)",
          "muEff": r"$\mu_{\mathrm{eff}}$ interpolation (SA)",
          "kappaEff": r"$\kappa_{\mathrm{eff}}$ interpolation (SA)"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--font-dir", type=Path, default=Path(
        "/mnt/c/texlive/2026/texmf-dist/fonts/opentype/public/tex-gyre"))
    args = parser.parse_args()
    family = "DejaVu Sans"
    if args.font_dir.is_dir():
        for face in ("regular", "italic", "bold", "bolditalic"):
            font_manager.fontManager.addfont(args.font_dir / f"texgyreheros-{face}.otf")
        family = "TeX Gyre Heros"
    plt.rcParams.update({
        "font.family": family, "font.size": 10, "axes.labelsize": 10,
        "axes.titlesize": 10, "xtick.labelsize": 10, "ytick.labelsize": 10,
        "legend.fontsize": 10, "mathtext.fontset": "cm", "pdf.fonttype": 3,
        "axes.linewidth": 0.65, "axes.spines.top": False,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })
    provenance = json.loads((HERE / "provenance.json").read_text())
    data = HERE / provenance["retained_data"]
    assert hashlib.sha256(data.read_bytes()).hexdigest() == provenance["retained_data_sha256"]
    limits = pd.read_csv(data, float_precision="round_trip")
    fig, axes = plt.subplots(1, 2, figsize=(170 / 25.4, 98 / 25.4),
                             constrained_layout=True)
    fig.get_layout_engine().set(w_pad=0.07, h_pad=0.08, wspace=0.04)
    for ax, fields, title in zip(
        axes, [("T", "U"), ("muEff", "kappaEff")],
        ["Gradient differencing", "Linear coefficient interpolation"],
    ):
        for field in fields:
            selected = limits[limits.component == field]
            q = selected.y_face_plus.to_numpy()
            h = selected.h_plus.to_numpy()
            line, = ax.plot(q, h, color=COLORS[field], label=LABELS[field], lw=1.2)
            np.testing.assert_array_equal(line.get_xdata(), q)
            np.testing.assert_array_equal(line.get_ydata(), h)
        ax.plot(q, q, color="#707D87", lw=0.8,
                label=r"Geometry cap: $h^+=y_f^+$")
        ax.set(xscale="log", yscale="log", xlim=provenance["xlim"],
               ylim=provenance["ylim"], title=title,
               xlabel=r"Face location $y_f^+$" + "\n(DNS wall units)",
               ylabel=r"Full centre separation $h^+$")
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_locator(ticker.LogLocator(base=10, numticks=12))
            axis.set_minor_locator(ticker.LogLocator(
                base=10, subs=np.arange(2, 10), numticks=100))
            axis.set_minor_formatter(ticker.NullFormatter())
        ax.grid(True, which="major", color="#D0D0D0", linewidth=0.5)
        ax.grid(True, which="minor", color="#ECECEC", linewidth=0.35)
        ax.legend(loc="upper left", frameon=False, handlelength=1.5,
                  handletextpad=0.5, borderaxespad=0.3, labelspacing=0.4)
    fig.canvas.draw()
    stem = HERE.parent / "M6Tw025_3pct_envelope"
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=300)
    plt.close(fig)
    print(f"Rendered {stem.name}: 170 x 98 mm, 10 pt, exact retained coordinates.")


if __name__ == "__main__":
    main()
