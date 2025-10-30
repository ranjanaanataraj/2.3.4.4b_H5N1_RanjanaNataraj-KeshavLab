#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
viz_kde3d.py — 3D ridge-style KDE for discrete groups.

Input CSV must contain:
  - value: numeric (e.g., GLS_count or Stalk_length)
  - group: categorical (clade/host/etc.)

Example:
  python viz_kde3d.py \
    --csv example_data/gls_demo.csv \
    --xcol GLS_count --gcol Group \
    --xmin 3 --xmax 11 --bw 0.4 --grid 400 \
    --elev 8 --azim -80 \
    --out bubble_kde3d.png
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def parse_args():
    p = argparse.ArgumentParser(description="3D KDE ridges over groups")
    p.add_argument("--csv", required=True, help="Input CSV (tidy)")
    p.add_argument("--xcol", default="value", help="Numeric column for KDE (default: value)")
    p.add_argument("--gcol", default="group", help="Group column (default: group)")
    p.add_argument("--xmin", type=float, default=None)
    p.add_argument("--xmax", type=float, default=None)
    p.add_argument("--grid", type=int, default=300, help="Grid points for KDE (default: 300)")
    p.add_argument("--bw", type=float, default=None,
                   help="Bandwidth factor for gaussian_kde (e.g., 0.4). If None, use Scott's rule.")
    p.add_argument("--alpha", type=float, default=0.6, help="Fill alpha (default: 0.6)")
    p.add_argument("--gap", type=float, default=1.0, help="Spacing between ridges along Y (default: 1.0)")
    p.add_argument("--elev", type=float, default=8.0, help="View elevation (default: 8)")
    p.add_argument("--azim", type=float, default=-80.0, help="View azimuth (default: -80)")
    p.add_argument("--out", required=True, help="Output image path (.png/.svg/.tiff)")
    return p.parse_args()

def kde_1d(x, grid, bw=None):
    kde = gaussian_kde(x, bw_method=bw)
    return kde(grid)

def shade_under_curve(ax, x, y_level, z, color, alpha):
    """Fill vertical polygon between KDE curve and z=0 plane at fixed y."""
    verts = [(x[0], y_level, 0.0)]
    verts += [(xi, y_level, zi) for xi, zi in zip(x, z)]
    verts += [(x[-1], y_level, 0.0)]
    poly = Poly3DCollection([verts], facecolor=color, alpha=alpha, linewidths=0)
    ax.add_collection3d(poly)

def main():
    a = parse_args()
    df = pd.read_csv(a.csv)
    if a.xcol not in df or a.gcol not in df:
        raise SystemExit(f"Columns {a.xcol!r} and/or {a.gcol!r} not found in {a.csv}")

    # Clean and bounds
    df = df[[a.xcol, a.gcol]].dropna()
    x = df[a.xcol].astype(float)
    xmin = a.xmin if a.xmin is not None else np.nanpercentile(x, 0.5)
    xmax = a.xmax if a.xmax is not None else np.nanpercentile(x, 99.5)
    grid = np.linspace(xmin, xmax, a.grid)

    groups = sorted(df[a.gcol].dropna().unique())
    n = len(groups)
    if n == 0:
        raise SystemExit("No groups found.")

    # Color cycle from Matplotlib
    colors = plt.rcParams["axes.prop_cycle"].by_key().get("color", ["C0","C1","C2","C3","C4"])

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_box_aspect([3, 1, 2])

    # Plot each group as a ridge along Y
    for i, g in enumerate(groups):
        sub = df[df[a.gcol] == g][a.xcol].values
        if len(sub) < 3:
            # too few to KDE; draw a tiny spike to indicate presence
            dens = np.zeros_like(grid)
        else:
            dens = kde_1d(sub, grid, bw=a.bw)

        y_level = i * a.gap
        color = colors[i % len(colors)]

        # Filled area + curve
        shade_under_curve(ax, grid, y_level, dens, color, a.alpha)
        ax.plot(grid, np.full_like(grid, y_level), dens, color=color, lw=2, alpha=1.0, label=str(g))

    # Aesthetics
    ax.view_init(elev=a.elev, azim=a.azim)
    ax.set_xlabel("Value", labelpad=10)
    ax.set_ylabel("")  # y is categorical groups
    ax.set_zlabel("Density", labelpad=10)
    ax.set_yticks([i * a.gap for i in range(n)])
    ax.set_yticklabels(groups)
    ax.xaxis._axinfo["grid"].update(color="gray", linestyle="dashed", linewidth=0.5)
    ax.zaxis._axinfo["grid"].update(color="gray", linestyle="dashed", linewidth=0.5)

    # Legend outside
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, title=a.gcol)

    plt.tight_layout()
    plt.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    main()
