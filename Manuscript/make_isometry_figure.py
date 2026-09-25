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
  (Jirinec et al. 2021); it is drawn qualitatively because that study reports
  species-level trends, not a community-level rate comparable with ours.

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

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 5.6), gridspec_kw={"width_ratios": [1.15, 1.0]})

# ---------------------------------------------------------------------------
# Panel A: rates and isometry contrast
# ---------------------------------------------------------------------------
rows = [
    ("wing", "Wing length (ΔL)", "fitted on log wing length", "#2B6CB0", "o", "#2B6CB0", "-"),
    ("mass", "Body mass (ΔM)", "fitted on log body mass", "#C53030", "o", "#C53030", "-"),
    ("expected_isometric_mass", "Mass change expected under isometry",
     "M ∝ L³: three times the wing slope", "#B7791F", "s", "white", "--"),
    ("contrast", "Isometry contrast (ln M − 3 ln L)", "fitted directly · full adjustment", "#2F855A", "D", "#2F855A", "-"),
    ("contrast_parsimonious", "Isometry contrast (ln M − 3 ln L)", "fitted directly · parsimonious adjustment",
     "#2F855A", "D", "white", "-"),
]
ys = [4.4, 3.4, 2.4, 1.1, 0.2]
ax1.axvline(0, color="#718096", lw=1.0, ls="--", zorder=1)
ax1.axhspan(1.85, 4.95, color="#F7FAFC", zorder=0)
ax1.axhspan(-0.4, 1.65, color="#F0FFF4", zorder=0)
XL = -13.5
for (key, lab, sub, col, mk, mfc, ls), y in zip(rows, ys):
    r = iso.loc[key]
    ax1.plot([r.lower, r.upper], [y, y], color=col, lw=1.8, ls=ls, zorder=3)
    ax1.plot(r.pct_per_decade, y, marker=mk, ms=7, color=col, markerfacecolor=mfc, markeredgewidth=1.5, zorder=4)
    ax1.text(XL, y + 0.13, lab, fontsize=8.2, fontweight="bold", color="#1A202C", va="center")
    ax1.text(XL, y - 0.19, sub, fontsize=7.2, color="#4A5568", va="center")
    ax1.text(r.upper + 0.35, y, f"{sgn(r.pct_per_decade)}%  [{sgn(r.lower)}, {sgn(r.upper)}]",
             fontsize=7.6, color=col, va="center", fontweight="bold")
    if key.startswith("contrast"):
        ax1.text(r.upper + 0.35, y - 0.28, f"P(β > 0) = {1 - r.p_neg:.2f}", fontsize=7, color=col, va="center")
ax1.text(XL, 4.85, f"Atlantic Forest, shared records (N = {N:,}; {len(spp)} species)", fontsize=8.5,
         fontweight="bold", color="#2D3748")
ax1.text(XL, 1.55, "Isometry contrast: > 0 = mass retained relative to wing", fontsize=8.5,
         fontweight="bold", color="#22543D")
ax1.set_xlim(XL - 0.4, 12.0); ax1.set_ylim(-0.5, 5.2)
ax1.set_yticks([]); ax1.spines["left"].set_visible(False)
ax1.set_xlabel("Rate of change (% per decade) · posterior mean and 95% CrI", fontsize=8.8)
ax1.set_title("A. Decadal rates and the isometry contrast", fontsize=10.5, fontweight="bold", loc="left", pad=10)

# ---------------------------------------------------------------------------
# Panel B: bivariate state space
# ---------------------------------------------------------------------------
pad = 1.0
xlim = (min(-6.0, spp.pct_per_decade_wing.min() - pad), max(4.0, spp.pct_per_decade_wing.max() + pad))
ylim = (min(-8.0, spp.pct_per_decade_mass.min() - pad), max(6.0, spp.pct_per_decade_mass.max() + pad))
ax2.set_xlim(xlim); ax2.set_ylim(ylim)
ax2.add_patch(Rectangle((0, ylim[0]), xlim[1], -ylim[0], facecolor="#FFF5F5", edgecolor="none", zorder=0))
ax2.text(xlim[1] - 0.15, ylim[0] + 0.4, "Longer wings, lower mass:\ndirection reported for the\ncentral Amazon (Jirinec et al. 2021)",
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
ax2.set_title("B. Bivariate state space", fontsize=10.5, fontweight="bold", loc="left", pad=10)
ax2.legend(handles=[
    Line2D([0], [0], color="#1A202C", lw=1.8, label="strict isometry (y = 3x)"),
    Line2D([0], [0], color="#718096", lw=1.1, ls=":", label="equal rates (y = x)"),
    Line2D([0], [0], marker="o", color="#276749", mec="white", ms=7, lw=1.8,
           label="assemblage estimate, marginal 95% CrIs"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor="#718096", markeredgecolor="#4A5568", alpha=0.6,
           ms=5, label=f"species posterior means (n = {len(spp)}; {n_decoupled} in the shaded quadrant)")],
    loc="upper left", frameon=True, facecolor="white", framealpha=0.9, edgecolor="#E2E8F0", fontsize=7.0)

fig.suptitle("Stan phylogenetic mixed models, fully adjusted (+ municipality)"
             + (f"; posterior pooled over {n_trees} trees" if n_trees else ""), fontsize=8.5, y=1.0)
plt.tight_layout()
fig_path = os.path.join(OUT, "fig-isometry.png")
fig.savefig(fig_path, bbox_inches="tight"); plt.close(fig)
print(f"Successfully created {fig_path} (N = {N}, {len(spp)} species, {n_decoupled} in the longer-wing/lower-mass quadrant)")
