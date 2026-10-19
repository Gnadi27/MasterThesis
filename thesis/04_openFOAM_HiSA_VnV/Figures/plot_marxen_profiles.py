"""Replot the Marxen downstream profiles with DejaVu Sans and 25% larger text.

The original constants, digitized reference arrays and normalization function
are read directly from ValData/DNS/Marxen/plot_marxen_profileinlet.py. Original
scripts and CSVs remain unchanged. Run from any directory; --data-root can
override the adjacent Files/ValData archive. Requires matplotlib/numpy/pandas.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.text import Text
import numpy as np


HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[2]
STEM = "04_Marxen_lamBL_profiles"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_reader(script):
    names = {"T_E", "U_E", "MU_E", "RHO_E", "NU_E", "ref_eta_t", "ref_t", "ref_eta_u", "ref_u"}
    tree = ast.parse(script.read_text(), filename=str(script))
    nodes = []
    for node in tree.body:
        if isinstance(node, ast.Import) and all(a.name in {"numpy", "pandas"} for a in node.names):
            nodes.append(node)
        elif isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id in names for t in node.targets):
            nodes.append(node)
        elif isinstance(node, ast.FunctionDef) and node.name == "load_and_normalize":
            nodes.append(node)
    namespace = {}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(script), "exec"), namespace)
    return namespace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=WORKSPACE.parents[1]/"Files/ValData")
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    source = args.data_root.resolve()/"DNS/Marxen/plot_marxen_profileinlet.py"
    scientific = source_reader(source)
    csvs = [source.parent/"profile_inlet/Run1"/f"x0p{x}.csv" for x in (25713, 32543, 40176, 48613)]
    stations = sorted([scientific["load_and_normalize"](p) for p in csvs], key=lambda s:s["x"])
    scale = 1.25
    plt.rcParams.update({
        "font.family":"DejaVu Sans", "mathtext.fontset":"dejavusans", "text.usetex":False,
        "font.size":12*scale, "axes.titlesize":14*scale, "axes.labelsize":13*scale,
        "xtick.labelsize":11*scale, "ytick.labelsize":11*scale,
        "legend.fontsize":8.5*scale, "axes.linewidth":.8,
        "figure.facecolor":"white", "savefig.facecolor":"white", "pdf.fonttype":42,
    })
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=180)
    for s, color in zip(stations, plt.get_cmap("viridis")(np.linspace(.15,.85,len(stations)))):
        label = rf"OpenFOAM, $x={s['x']:.3f}$ m, $\sqrt{{Re_x}}={s['sqrt_Re_x']:.0f}$"
        for ax, key in zip(axes, ("U_star", "T_star")):
            ax.plot(s[key], s["Eta"], lw=2, color=color, label=label)
    style = dict(color="black", lw=1, marker="o", ms=5, mfc="white", label="Marxen et al. Fig. 11")
    axes[0].plot(scientific["ref_u"], scientific["ref_eta_u"], **style)
    axes[1].plot(scientific["ref_t"], scientific["ref_eta_t"], **style)
    for i, ax in enumerate(axes):
        ax.set_xlabel(r"$U/u_e$" if i == 0 else r"$T/T_e$")
        ax.set_ylabel(r"$\eta=y/L_r(x),\quad L_r=\sqrt{\nu_e\,x/u_e}$" if i == 0 else r"$\eta=y/L_r(x)$")
        ax.set_title("(a) Velocity profile" if i == 0 else "(b) Temperature profile")
        ax.set_xlim((0,1.05) if i == 0 else (.9,6.5))
        ax.set_ylim(0,25)
        ax.set_xticks(np.arange(0,1.01,.2) if i == 0 else np.arange(1,7))
        ax.set_yticks(np.arange(0,26,5))
        ax.grid(True, alpha=.3)
        ax.legend(loc="upper left" if i == 0 else "upper right", framealpha=.95)
    fig.tight_layout(pad=.7, w_pad=.8)
    lines = []
    for ax in axes:
        for line in ax.lines:
            xy = np.column_stack((line.get_xdata(), line.get_ydata())).astype("<f8")
            lines.append({"label":line.get_label(), "points":len(xy),
                          "sha256":hashlib.sha256(xy.tobytes()).hexdigest()})
    audit = WORKSPACE/"tmp/laminar_font_review/render_audit.json"
    checked = False
    if audit.exists():
        expected = json.loads(audit.read_text())["marxen_downstream"]["exports"][0]["lines"]
        assert [p["sha256"] for p in lines] == [p["sha256"] for p in expected], "Plotted values changed"
        checked = True
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for text in fig.findobj(Text):
        if text.get_visible() and text.get_text():
            box = text.get_window_extent(renderer)
            assert box.x0 >= -.5 and box.y0 >= -.5 and box.x1 <= fig.bbox.width+.5 and box.y1 <= fig.bbox.height+.5, text.get_text()
    for ext in ("png", "pdf"):
        fig.savefig(args.output_dir/f"{STEM}.{ext}", dpi=180)
    plt.close(fig)
    provenance = {"font_family":"DejaVu Sans", "font_scale":scale,
                  "source_script":str(source), "source_script_sha256":sha(source),
                  "sources":[{"path":str(p),"sha256":sha(p)} for p in csvs],
                  "lines":lines, "matches_previous_render_audit":checked,
                  "normalization":"original source load_and_normalize, unchanged",
                  "text_inside_figure":True}
    (args.output_dir/"marxen_profile_font_provenance.json").write_text(json.dumps(provenance,indent=2)+"\n")
    print(f"Saved {STEM}; historical plotted-value verification: {checked}")


if __name__ == "__main__":
    main()
