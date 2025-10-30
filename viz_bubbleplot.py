#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
viz_bubbleplot.py — Bubble plot for H5 GLS × NA subtype.

Input CSV must contain at least:
    #GLS          → numeric glycosylation counts
    NA            → NA subtype string (e.g. "N1")
    Stalk_length  → numeric stalk length (for color scale)
    Frequency     → numeric (bubble area)
"""

import argparse, numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns

def parse_args():
    p = argparse.ArgumentParser(description="Bubble plot: H5 GLS × NA subtype")
    p.add_argument("--csv", required=True, help="Input CSV with #GLS, NA, Stalk_length, Frequency")
    p.add_argument("--out", required=True, help="Output image (pdf/png)")
    p.add_argument("--title", default="H5 GLS × NA stalk length — counts")
    return p.parse_args()

def size_scale(freq):
    fmin, fmax = float(freq.min()), float(freq.max())
    s_min, s_max = 10.0, 2500.0
    def fn(v):
        if fmax == fmin:
            return (s_min + s_max) / 2.0
        return s_min + (np.sqrt(v - fmin) / np.sqrt(fmax - fmin)) * (s_max - s_min)
    return np.array([fn(v) for v in freq])

def main():
    a = parse_args()
    df = pd.read_csv(a.csv)

    # --- clean inputs ---
    df["#GLS"] = pd.to_numeric(df["#GLS"], errors="coerce")
    df["NA_norm"] = (df["NA"].astype(str).str.upper()
                     .str.replace(r"\s+", "", regex=True)
                     .str.extract(r"(N\d+)", expand=False))
    df = df.dropna(subset=["#GLS", "NA_norm", "Frequency", "Stalk_length"])

    gls_vals = sorted(df["#GLS"].dropna().unique().astype(int))
    xmap = {v: i for i, v in enumerate(gls_vals)}
    x = df["#GLS"].map(xmap)

    na_order = [f"N{i}" for i in range(1, 10)]
    present = [n for n in na_order if n in df["NA_norm"].unique()]
    ymap = {n: i for i, n in enumerate(present)}
    y = df["NA_norm"].map(ymap)

    sizes = size_scale(df["Frequency"].astype(float).to_numpy())

    plt.rcdefaults()
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "font.family": "sans-serif",
        "font.size": 12,
        "axes.labelsize": 13,
        "axes.titlesize": 14,
    })

    fig, ax = plt.subplots(figsize=(12, 8))
    sc = ax.scatter(x, y, s=sizes, c=df["Stalk_length"], cmap="viridis",
                    alpha=0.9, edgecolors="black", linewidths=0.6, zorder=3)

    # styling
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_color("#B0B0B0"); ax.spines[s].set_linewidth(0.8)
    ax.grid(True, axis="y", color="#E6E6E6", linewidth=0.8)
    ax.set_axisbelow(True)

    ax.set_xticks(range(len(gls_vals)))
    ax.set_xticklabels([str(v) for v in gls_vals])
    ax.set_xlabel("H5 #GLS")
    ax.set_yticks(range(len(present)))
    ax.set_yticklabels(present)
    ax.set_ylabel("NA subtype")
    ax.set_title(a.title, pad=10)

    cbar = plt.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("Stalk length (aa)")

    freq = df["Frequency"].astype(float).to_numpy()
    ref = np.unique(np.percentile(freq, [25, 50, 90]).astype(int))
    handles = [plt.scatter([], [], s=size_scale(np.array([v]))[0],
                           facecolors="none", edgecolors="black", linewidths=0.8)
               for v in ref]
    ax.legend(handles, [str(v) for v in ref], title="Frequency",
              frameon=False, scatterpoints=1,
              bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.tick_params(length=0)

    plt.tight_layout()
    plt.savefig(a.out, bbox_inches="tight")
    print(f"[bubbleplot] Saved {a.out}")

if __name__ == "__main__":
    main()
