#!/usr/bin/env python3
"""
make_isometry_figure.py
-----------------------
Wing vs body mass temporal allometry and the geometric isometry test on shared
records (Manuscript/images/fig-isometry.png), drawn from the Stan phylogenetic
tier (posterior pooled over trees), the engine every estimate in the text uses.

Panel A: decadal rates (% per decade, posterior mean and 95% CrI) under the
  fully adjusted set (+ municipality): wing length, body mass, the mass change
  expected under strict isometry (M ~ L^3, i.e. three times the log wing slope),
  and the directly fitted isometry contrast ln M - 3 ln L (full and parsimonious
  adjustment sets).
Panel B: bivariate state space (dWing vs dMass). Wing and mass are fitted in
  separate models, so there is no joint posterior: the assemblage estimate is
  shown with its two marginal 95% CrIs, not an ellipse. Species points are
  posterior means of fixed + species deviation. The quadrant of longer wings
  and lower mass is the direction of change reported for the central Amazon
  (Jirinec et al. 2021), drawn as a direction; the rates themselves are compared
  in fig-amazon (Analysis/scripts/jirinec_benchmark.R).

Inputs (Rscript Analysis/scripts/make_figures.R export):
  Analysis/output/figure_data/fig_isometry_stan.csv
  Analysis/output/figure_data/fig_isometry_species_stan.csv
  Analysis/output/figure_data/fig_scalars.csv            (stan_n_trees)
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "images")
FD = os.path.normpath(os.path.join(HERE, "..", "Analysis", "output", "figure_data"))

iso = pd.read_csv(os.path.join(FD, "fig_isometry_stan.csv")).set_index("quantity")
spp = pd.read_csv(os.path.join(FD, "fig_isometry_species_stan.csv"))
scal = pd.read_csv(os.path.join(FD, "fig_scalars.csv")).set_index("key")["value"]
n_trees = int(float(scal["stan_n_trees"])) if "stan_n_trees" in scal else None
N = int(iso.loc["contrast", "N"])

plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 300})

def neg(val, fmt=".1f"):
    return f"{val:{fmt}}".replace("-", "−")

def sgn(val, fmt=".1f"):
    return neg(val, "+" + fmt)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.6), gridspec_kw={"width_ratios": [1.0, 1.0]})

# ---------------------------------------------------------------------------
# Panel A: rates and isometry contrast -- a plain forest plot (Shinichi Nakagawa's
# comment, 2026-09-30: the earlier version with shaded blocks, sub-labels and value
# annotations was hard to read). Values are given in the text and legend.
# ---------------------------------------------------------------------------
rows = [
    ("wing", "Wing length", "#4A5568", "o", "#4A5568", "-"),
    ("mass", "Body mass", "#4A5568", "o", "#4A5568", "-"),
    ("expected_isometric_mass", "Mass expected under isometry", "#4A5568", "s", "white", "--"),
    ("contrast", "Isometry contrast, full adjustment", "#2F855A", "D", "#2F855A", "-"),
    ("contrast_parsimonious", "Isometry contrast, parsimonious", "#2F855A", "D", "white", "-"),
]
ys = [4, 3, 2, 0.8, -0.2]
ax1.axvline(0, color="#A0AEC0", lw=0.9, ls="--", zorder=1)
ax1.axhline(1.4, color="#E2E8F0", lw=0.8, zorder=0)
for (key, lab, col, mk, mfc, ls), y in zip(rows, ys):
    r = iso.loc[key]
    ax1.plot([r.lower, r.upper], [y, y], color=col, lw=1.8, ls=ls, zorder=3)
    ax1.plot(r.pct_per_decade, y, marker=mk, ms=7, color=col, markerfacecolor=mfc, markeredgewidth=1.5, zorder=4)
ax1.set_yticks(ys); ax1.set_yticklabels([r[1] for r in rows], fontsize=8.6)
ax1.tick_params(axis="y", length=0)
lo = min(iso.loc[[r[0] for r in rows], "lower"].min(), -1.0); hi = max(iso.loc[[r[0] for r in rows], "upper"].max(), 1.0)
ax1.set_xlim(lo - 0.8, hi + 0.8); ax1.set_ylim(-0.8, 4.6)
ax1.spines["left"].set_visible(False)
ax1.set_xlabel("Change (% per decade)", fontsize=8.8)
ax1.set_title("A", fontsize=10.5, fontweight="bold", loc="left", pad=10)

# ---------------------------------------------------------------------------
# Panel B: bivariate state space
# ---------------------------------------------------------------------------
pad = 1.0
xlim = (min(-6.0, spp.pct_per_decade_wing.min() - pad), max(4.0, spp.pct_per_decade_wing.max() + pad))
ylim = (min(-8.0, spp.pct_per_decade_mass.min() - pad), max(6.0, spp.pct_per_decade_mass.max() + pad))
ax2.set_xlim(xlim); ax2.set_ylim(ylim)
ax2.add_patch(Rectangle((0, ylim[0]), xlim[1], -ylim[0], facecolor="#FFF5F5", edgecolor="none", zorder=0))
ax2.text(xlim[1] - 0.15, ylim[0] + 0.4, "longer wings, lower mass\n(central Amazon)",
         fontsize=7.4, color="#9B2C2C", style="italic", ha="right", va="bottom")
ax2.axvline(0, color="#A0AEC0", lw=0.9, ls="--", zorder=1)
ax2.axhline(0, color="#A0AEC0", lw=0.9, ls="--", zorder=1)
xv = np.linspace(xlim[0], xlim[1], 100)
ax2.plot(xv, 3 * xv, color="#1A202C", lw=1.8, zorder=2)
ax2.plot(xv, xv, color="#718096", lw=1.1, ls=":", zorder=2)
ax2.scatter(spp.pct_per_decade_wing, spp.pct_per_decade_mass, s=24, color="#718096", alpha=0.45,
            edgecolors="#4A5568", linewidths=0.5, zorder=3)
w, m = iso.loc["wing"], iso.loc["mass"]
ax2.errorbar(w.pct_per_decade, m.pct_per_decade,
             xerr=[[w.pct_per_decade - w.lower], [w.upper - w.pct_per_decade]],
             yerr=[[m.pct_per_decade - m.lower], [m.upper - m.pct_per_decade]],
             fmt="o", ms=8, color="#276749", mec="white", mew=1.3, ecolor="#276749", elinewidth=1.8,
             capsize=3.5, zorder=5)
n_decoupled = int(((spp.pct_per_decade_wing > 0) & (spp.pct_per_decade_mass < 0)).sum())
ax2.set_xlabel("Wing length change (% per decade)", fontsize=8.8)
ax2.set_ylabel("Body mass change (% per decade)", fontsize=8.8)
ax2.set_title("B", fontsize=10.5, fontweight="bold", loc="left", pad=10)
ax2.legend(handles=[
    Line2D([0], [0], color="#1A202C", lw=1.8, label="isometry (y = 3x)"),
    Line2D([0], [0], color="#718096", lw=1.1, ls=":", label="equal rates (y = x)"),
    Line2D([0], [0], marker="o", color="#276749", mec="white", ms=7, lw=1.8,
           label="assemblage, 95% CrIs"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor="#718096", markeredgecolor="#4A5568", alpha=0.6,
           ms=5, label=f"species (n = {len(spp)})")],
    loc="upper left", frameon=True, facecolor="white", framealpha=0.9, edgecolor="#E2E8F0", fontsize=7.4)

plt.tight_layout()
fig_path = os.path.join(OUT, "fig-isometry.png")
fig.savefig(fig_path, bbox_inches="tight"); plt.close(fig)
print(f"Successfully created {fig_path} (N = {N}, {len(spp)} species, {n_decoupled} in the longer-wing/lower-mass quadrant)")
