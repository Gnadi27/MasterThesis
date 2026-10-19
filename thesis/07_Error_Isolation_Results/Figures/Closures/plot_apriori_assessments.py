#!/usr/bin/env python3
"""Replot the twelve thesis a-priori assessments from saved solver sweeps.

Requires Python 3, matplotlib, numpy and pandas. Run without arguments to
use the byte-for-byte CSV copies in apriori_sources/. The manifest records
their original WSL paths, hashes and analysis scripts. No closure is solved
again, and the saved error columns and validity filters are unchanged.

Each output has shear on the left and heat flux on the right. The fixed
8-inch canvas and original 10.5-point font preserve the printed lettering
at the manuscript's existing inclusion widths, while the shorter canvas
allows a Li/Griffin pair and both captions to share a page.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "thesis-apriori-mpl"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
SOURCES = HERE / "apriori_sources"


def read_sweep(source: dict) -> pd.DataFrame:
    path = SOURCES / source["copy"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != source["sha256"]:
        raise ValueError(f"Saved sweep has changed: {path}")
    data = pd.read_csv(path, float_precision="round_trip")
    validity = data[source["validity_column"]]
    if not pd.api.types.is_bool_dtype(validity):
        raise ValueError(f"Expected boolean validity column in {path}")
    return data.loc[validity].sort_values(source["height_column"])


def draw_figure(spec: dict) -> dict:
    sweeps = [(source, read_sweep(source)) for source in spec["sources"]]
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.1))
    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.22, top=0.85, wspace=0.32)
    columns = ("tau_w_error_percent", "q_w_error_percent")
    colors = ("#1769aa", "#b42318")
    labels = (r"$\tau_w$ error [%]", r"$q_w$ error [%]")
    heights = np.concatenate([data[source["height_column"]] for source, data in sweeps])
    audit = {"output": spec["output"], "panels": []}
    for ax, column, color, label in zip(axes, columns, colors, labels):
        panel = {"column": column, "series": []}
        ax.axhline(0.0, color="0.3", linewidth=1.0)
        for source, data in sweeps:
            finite = data.loc[np.isfinite(data[column])]
            edge = source["role"] == "edge"
            if not finite.empty:
                line, = ax.plot(
                    finite[source["height_column"]], finite[column], color=color,
                    marker=None if edge else "o", linestyle="--" if edge else "-",
                    alpha=0.35 if edge else 1.0, linewidth=1.8 if edge else 2.0,
                    zorder=2 if edge else 3,
                )
                # Verify the plotted points, including the edge-only endpoint,
                # directly against each source's unmodified saved columns.
                np.testing.assert_array_equal(line.get_xdata(), finite[source["height_column"]])
                np.testing.assert_array_equal(line.get_ydata(), finite[column])
            panel["series"].append({
                "role": source["role"], "point_count": len(finite),
                "x": finite[source["height_column"]].tolist(),
                "y": finite[column].tolist(),
            })
        if not any(series["point_count"] for series in panel["series"]):
            # Zero reference heat flux has no percentage error. Do not plot
            # zeros as if they were valid heat-flux predictions.
            message = "Reference $q_w=0$\nPercentage error undefined"
            if spec["model"] == "Griffin":
                data = sweeps[0][1]
                low, high = data["q_w_model_W_per_m2"].agg(["min", "max"])
                message += f"\nModeled $q_w$: {low:.4g} to {high:.4g} W/m$^2$"
            ax.text(0.5, 0.62, message, ha="center", va="center", transform=ax.transAxes)
        ax.set_xscale("log")
        ax.set_xlim(0.9 * heights.min(), 1.1 * heights.max())
        ax.set_xlabel(r"Hypothetical first-cell $y_1^+$" if spec["model"] == "Li"
                      else r"Matching location $y_m^+$")
        ax.set_ylabel(label)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        audit["panels"].append(panel)
    if spec["model"] == "Li":
        handles = [
            Line2D([], [], color="0.25", marker="o", linewidth=2.0, label="Two-point"),
            Line2D([], [], color="0.25", linestyle="--", alpha=0.35, linewidth=1.8,
                   label="Boundary-layer edge"),
        ]
        fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.55, 1.0),
                   ncol=2, frameon=False)
    output = HERE.parent / spec["output"]
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in axes:
        bounds = ax.get_tightbbox(renderer)
        if (bounds.x0 < 0 or bounds.y0 < 0 or bounds.x1 > fig.bbox.width
                or bounds.y1 > fig.bbox.height):
            raise ValueError(f"Text extends beyond canvas: {output}")
    for suffix in (".png", ".pdf"):
        fig.savefig(output.with_suffix(suffix), dpi=220)
    plt.close(fig)
    return audit


def main() -> None:
    plt.rcParams.update({
        "font.size": 10.5, "axes.labelsize": 10.5, "legend.fontsize": 9.5,
        "axes.grid": True, "grid.alpha": 0.24, "grid.linewidth": 0.7,
        "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 2.0, "lines.markersize": 6.0,
        "pdf.fonttype": 42, "savefig.bbox": None,
    })
    manifest = json.loads((SOURCES / "manifest.json").read_text())
    audit = [draw_figure(spec) for spec in manifest["figures"]]
    (SOURCES / "plot_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(f"Replotted {len(audit)} figures; all plotted values match the saved sweeps.")


if __name__ == "__main__":
    main()
