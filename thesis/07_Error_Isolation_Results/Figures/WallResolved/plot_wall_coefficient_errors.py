#!/usr/bin/env python3
r"""Replot Figures 7.2/7.3 at their printed size from retained numerical data.

Run with Python, numpy, pandas and matplotlib. The adjacent CSV retains the
curated values used in the previous figures; the JSON records WSL provenance.
The PDF canvas is exactly the thesis text width (170 mm), so 9.5-point labels
remain 9.5-point labels when included at width=\linewidth. PNGs are also saved.
"""
from pathlib import Path
import argparse
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, ticker
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
STYLES = {
    "Martin-Duan Part 2": ("^", "#1f77b4", "Martin–Duan Part 2"),
    "Martin-Duan Part 4": ("s", "#d95f02", "Martin–Duan Part 4"),
    "Zhang": ("o", "#2ca02c", "Zhang"),
    "Huang": ("*", "#9467bd", "Huang (M11Tw020)"),
}
PANELS = (
    ("minus_Bq_DNS", "Cf", "(a) Skin-friction error", (-0.008, 0.202)),
    ("minus_Bq_DNS", "Ch", "(b) Heat-transfer error", (-0.008, 0.202)),
    ("M_tau_DNS", "Cf", "(c) Skin-friction error", (0.074, 0.198)),
    ("M_tau_DNS", "Ch", "(d) Heat-transfer error", (0.074, 0.198)),
)


def configure_fonts(font_dir):
    # Heros is the OpenType equivalent of the thesis's Helvetica/Nimbus Sans.
    if font_dir.is_dir():
        for face in ("regular", "italic", "bold", "bolditalic"):
            font_manager.fontManager.addfont(font_dir / f"texgyreheros-{face}.otf")
        family = "TeX Gyre Heros"
    else:
        family = "DejaVu Sans"
    plt.rcParams.update({
        "font.family": family, "font.size": 9.5,
        "axes.labelsize": 9.5, "axes.titlesize": 10,
        "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
        "legend.fontsize": 9.5, "mathtext.fontset": "cm",
        "axes.linewidth": 0.6, "pdf.fonttype": 3,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })


def load_data():
    path = HERE / "wall_coefficient_error_data.csv"
    provenance = json.loads((HERE / "wall_coefficient_error_provenance.json").read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == provenance["data_sha256"]
    data = pd.read_csv(path, float_precision="round_trip")
    assert len(data) == 28
    for coefficient in ("Cf", "Ch"):
        rows = data.loc[data[f"{coefficient}_status"].eq("Compared")]
        expected = 100 * (rows[f"{coefficient}_RANS"] / rows[f"{coefficient}_DNS"] - 1)
        np.testing.assert_allclose(rows[f"{coefficient}_error_pct"], expected,
                                   rtol=1e-11, atol=1e-11)
    for column in ("minus_Bq_DNS", "M_tau_DNS"):
        pair = data.pivot(index="case_label", columns="model", values=column)
        np.testing.assert_allclose(pair.SA, pair.kwSST, rtol=0, atol=1e-12)
    return data


def annotate_cases(ax, rows, xcolumn, coefficient):
    """Place full case labels without covering symbols or other labels.

    All offsets are in display points and only affect annotations. A fine
    leader line connects displaced labels to their exact scatter coordinates.
    """
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    scale = fig.dpi / 72
    points = [(float(getattr(row, xcolumn)), float(getattr(row, f"{coefficient}_error_pct")))
              for row in rows.itertuples(index=False)]
    display = ax.transData.transform(points)
    point_boxes = [Bbox.from_bounds(x - 4*scale, y - 4*scale, 8*scale, 8*scale)
                   for x, y in display]
    occupied = []
    labels = []
    # Long labels first; deterministic ties retain source/case order.
    order = sorted(range(len(rows)), key=lambda i: -len(rows.iloc[i].case_label))
    for i in order:
        case = rows.iloc[i].case_label
        annotation = ax.annotate(case, points[i], xytext=(6, 5),
                                 textcoords="offset points", va="center", fontsize=9.5,
                                 color="#222222", zorder=5)
        candidates = []
        for dy in (7, -7, 0, 14, -14, 21, -21, 28, -28):
            for dx, ha in ((6, "left"), (-6, "right"), (0, "center"),
                           (14, "left"), (-14, "right"), (24, "left"), (-24, "right")):
                annotation.set_position((dx, dy))
                annotation.set_ha(ha)
                box = annotation.get_window_extent(renderer).expanded(1.04, 1.15)
                outside = max(0, ax.bbox.x0-box.x0) + max(0, box.x1-ax.bbox.x1)
                outside += max(0, ax.bbox.y0-box.y0) + max(0, box.y1-ax.bbox.y1)
                overlaps = sum(box.overlaps(other) for other in occupied + point_boxes)
                score = outside*10000 + overlaps*100000 + np.hypot(dx, dy)
                candidates.append((score, dx, dy, ha, box))
        score, dx, dy, ha, box = min(candidates, key=lambda item: item[0])
        if score >= 100000:
            raise ValueError(f"Cannot place {case} without overlap in {coefficient}/{xcolumn}")
        annotation.set_position((dx, dy))
        annotation.set_ha(ha)
        occupied.append(box)
        if np.hypot(dx, dy) > 12:
            ax.annotate("", points[i], xytext=(dx, dy), textcoords="offset points",
                        arrowprops={"arrowstyle": "-", "color": "0.55", "lw": 0.45,
                                    "shrinkA": 1, "shrinkB": 4}, zorder=2)
        labels.append(annotation)
    return labels


def make_figure(data, model, output_dir):
    selected = data.loc[data.model.eq(model)]
    assert len(selected) == 14
    assert set(selected.loc[selected.Ch_status.ne("Compared"), "case_label"]) == {
        "M5T5", "M2p5", "M6Tw076"}
    fig, axes = plt.subplots(4, 1, figsize=(170 / 25.4, 208 / 25.4), dpi=150)
    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.060, top=0.925, hspace=0.50)
    handles = [Line2D([], [], linestyle="none", marker=m, color=c,
                      markeredgecolor="#222222", markeredgewidth=0.5,
                      markersize=7 if m == "*" else 5.5, label=label)
               for m, c, label in STYLES.values()]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.53, 0.998),
               ncol=2, frameon=False, columnspacing=1.8, handletextpad=0.5)
    audit = []
    for ax, (xcolumn, coefficient, title, limits) in zip(axes, PANELS):
        rows = selected.loc[selected[f"{coefficient}_status"].eq("Compared")]
        expected_count = 14 if coefficient == "Cf" else 11
        assert len(rows) == expected_count
        for row in rows.itertuples(index=False):
            marker, color, _ = STYLES[row.source]
            point = [[getattr(row, xcolumn), getattr(row, f"{coefficient}_error_pct")]]
            scatter = ax.scatter(*point[0], marker=marker, color=color,
                                 s=55 if marker == "*" else 30, edgecolor="#222222",
                                 linewidth=0.5, zorder=4)
            np.testing.assert_array_equal(scatter.get_offsets(), point)
        # Match the previous ordinate ranges, with extra room for large labels.
        values = rows[f"{coefficient}_error_pct"]
        pad = max(1.0, .03 * (values.max() - values.min()))
        ax.set_ylim(5*np.floor((values.min()-pad)/5), 5*np.ceil((values.max()+pad)/5))
        ax.set_xlim(limits)
        ax.set_title(title, loc="left", pad=4)
        ax.set_ylabel(rf"$\varepsilon_{coefficient[-1]}$ [%]", labelpad=4)
        if xcolumn == "minus_Bq_DNS":
            ax.set_xlabel(r"DNS wall heat-transfer rate, $-B_q$", labelpad=2)
            ax.set_xticks([0, .05, .10, .15, .20])
        else:
            ax.set_xlabel(r"DNS friction Mach number, $M_\tau=u_\tau/a_w$", labelpad=2)
            ax.set_xticks([.08, .10, .12, .14, .16, .18])
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
        ax.tick_params(length=3, width=.6, pad=2)
        ax.axhline(0, color="0.45", linewidth=.7, zorder=1)
        ax.grid(color="0.87", linewidth=.5)
        ax.set_axisbelow(True)
        annotate_cases(ax, rows, xcolumn, coefficient)
        audit.append({"x": xcolumn, "coefficient": coefficient,
                      "verified_points": len(rows), "cases": rows.case_label.tolist()})
    fig.canvas.draw()
    for ax in axes:
        box = ax.get_tightbbox(fig.canvas.get_renderer())
        assert fig.bbox.contains(box.x0, box.y0) and fig.bbox.contains(box.x1, box.y1)
    stem = output_dir / f"DNS_Cf_Ch_relative_error_{model}"
    fig.savefig(stem.with_suffix(".pdf"), metadata={"Creator": Path(__file__).name,
                                                 "Title": f"{model}: DNS-relative wall-coefficient errors"})
    fig.savefig(stem.with_suffix(".png"), dpi=300)
    plt.close(fig)
    return {"panels": audit, "width_mm": 170, "height_mm": 208,
            "label_font_pt": 9.5, "title_font_pt": 10}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    parser.add_argument("--font-dir", type=Path, default=Path(
        "/mnt/c/texlive/2026/texmf-dist/fonts/opentype/public/tex-gyre"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    configure_fonts(args.font_dir)
    data = load_data()
    audit = {model: make_figure(data, model, args.output_dir) for model in ("kwSST", "SA")}
    (args.output_dir / "wall_coefficient_layout_audit.json").write_text(json.dumps(audit, indent=2)+"\n")
    print("Replotted both figures; verified all 100 scatter points and full case labels.")


if __name__ == "__main__":
    main()
