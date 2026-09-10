#!/usr/bin/env python3
"""
make_isometry_figure.py
-----------------------
Generates the bivariate wing vs body mass temporal allometry and geometric isometry
figure for the Atlantic Forest passerine manuscript (saved to Manuscript/images/fig-isometry.png).

Panel A: Rate Estimates and Isometry Contrasts (% per decade)
  - M1a (Observer-only; parsimonious model avoiding 85.9% spatial redundancy with municipality):
    * Observed wing length rate (ΔL)
    * Observed body mass rate (ΔM)
    * Expected isometric mass rate (M ∝ L³ -> 3 * ΔL)
    * Formal isometry contrast (ΔM - 3ΔL): +4.85%/dec [0.81, 8.89%], z = 2.35, p = 0.019 (REJECTS ISOMETRY)
  - M1 (Both observer + municipality crossed; variance competition widens interval):
    * Isometry contrast: +4.07%/dec [-0.02, 8.16%], z = 1.95, p = 0.051 (borderline)
  - Central Amazon benchmark (Jirinec et al. 2021) showing morphological decoupling

Panel B: Bivariate Morphological State Space (ΔWing vs ΔMass)
  - Strict geometric isometry line y = 3x (M ∝ L³)
  - 95% bivariate confidence ellipses for M1a (primary) and M1 (comparison)
  - 72 individual species estimates (M3 controls)
  - Morphological decoupling quadrant (Amazonian pattern) vs Atlantic Forest mass retention
"""

import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle, Polygon
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "images")
os.makedirs(OUT, exist_ok=True)

DATA_DIR = os.path.normpath(os.path.join(HERE, "..", "Analysis", "output", "figure_data"))
SUMMARY_FILE = os.path.join(DATA_DIR, "fig_isometry_summary.json")
SPECIES_FILE = os.path.join(DATA_DIR, "fig_isometry_species.csv")

with open(SUMMARY_FILE, "r") as f:
    s = json.load(f)

spp_df = pd.read_csv(SPECIES_FILE)

# Styling
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
})

def neg(val, fmt=".2f"):
    txt = f"{val:{fmt}}"
    return txt.replace("-", "−")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 5.8), gridspec_kw={"width_ratios": [1.2, 1.25]})

# ===========================================================================
# PANEL A: Forest Plot of Rates and Isometry Contrasts
# ===========================================================================

items = [
    # M1a Parsimonious Tier
    {
        "y": 6.2,
        "label": "M1a: Observed wing length (ΔL)",
        "sublabel": "Observer control · 7,577 shared records",
        "est": s["m1a_wing_pct_decade"],
        "lo": s["m1a_wing_lower"],
        "hi": s["m1a_wing_upper"],
        "color": "#2B6CB0",
        "marker": "o",
        "ms": 7,
        "ls": "-",
        "mfc": "#2B6CB0"
    },
    {
        "y": 5.2,
        "label": "M1a: Observed body mass (ΔM)",
        "sublabel": "Observer control · stable body mass",
        "est": s["m1a_mass_pct_decade"],
        "lo": s["m1a_mass_lower"],
        "hi": s["m1a_mass_upper"],
        "color": "#C53030",
        "marker": "o",
        "ms": 7,
        "ls": "-",
        "mfc": "#C53030"
    },
    {
        "y": 4.2,
        "label": "M1a: Expected isometric mass",
        "sublabel": "M ∝ L³ (expected: 3 × ΔL)",
        "est": s["m1a_isometric_mass_expected"],
        "lo": s["m1a_isometric_mass_expected_lower"],
        "hi": s["m1a_isometric_mass_expected_upper"],
        "color": "#D69E2E",
        "marker": "s",
        "ms": 6.5,
        "ls": "--",
        "mfc": "white"
    },
    {
        "y": 3.0,
        "label": "M1a Isometry contrast (ΔM − 3ΔL)",
        "sublabel": "Parsimonious (z = 2.35, p = 0.019; excludes 0)",
        "est": s["m1a_iso_contrast"],
        "lo": s["m1a_iso_lower"],
        "hi": s["m1a_iso_upper"],
        "color": "#2F855A", # Green
        "marker": "D",
        "ms": 7.5,
        "ls": "-",
        "mfc": "#2F855A"
    },
    # M1 Comparison Tier
    {
        "y": 1.9,
        "label": "M1 Isometry contrast (ΔM − 3ΔL)",
        "sublabel": "Both controls (z = 1.95, p = 0.051; 86% site overlap)",
        "est": s["m1_iso_contrast"],
        "lo": s["m1_iso_lower"],
        "hi": s["m1_iso_upper"],
        "color": "#805AD5", # Purple
        "marker": "D",
        "ms": 6.5,
        "ls": ":",
        "mfc": "white"
    },
    # Amazon Benchmark
    {
        "y": 0.8,
        "label": "Amazon: Wing length rate",
        "sublabel": "Jirinec et al. (2021) · longer wings",
        "est": s["amazon_wing_pct_decade"],
        "lo": s["amazon_wing_lower"],
        "hi": s["amazon_wing_upper"],
        "color": "#4A5568",
        "marker": "^",
        "ms": 6.5,
        "ls": "-",
        "mfc": "#4A5568"
    },
    {
        "y": 0.0,
        "label": "Amazon: Body mass rate",
        "sublabel": "Jirinec et al. (2021) · systematic mass loss",
        "est": s["amazon_mass_pct_decade"],
        "lo": s["amazon_mass_lower"],
        "hi": s["amazon_mass_upper"],
        "color": "#9B2C2C",
        "marker": "v",
        "ms": 6.5,
        "ls": "-",
        "mfc": "#9B2C2C"
    }
]

# Reference zero line
ax1.axvline(0, color="#718096", lw=1.0, ls="--", zorder=1)

# Shaded grouping bands
ax1.axhspan(3.7, 6.7, color="#F7FAFC", alpha=0.9, zorder=0) # M1a data
ax1.axhspan(1.4, 3.5, color="#F0FFF4", alpha=0.9, zorder=0) # Isometry contrasts
ax1.axhspan(-0.6, 1.3, color="#FFF5F5", alpha=0.9, zorder=0) # Amazon

# Group headers
ax1.text(-13.8, 6.65, "M1a Parsimonious model (Observer only):", fontsize=8.5, fontweight="bold", color="#2D3748")
ax1.text(-13.8, 3.45, "Geometric Isometry Contrasts (ΔM − 3ΔL):", fontsize=8.5, fontweight="bold", color="#22543D")
ax1.text(-13.8, 1.25, "Central Amazon benchmark (Jirinec et al. 2021):", fontsize=8.5, fontweight="bold", color="#742A2A")

for it in items:
    y = it["y"]
    ax1.plot([it["lo"], it["hi"]], [y, y], color=it["color"], lw=1.8, ls=it["ls"], zorder=3)
    ax1.plot(it["est"], y, marker=it["marker"], markersize=it["ms"], color=it["color"],
             markerfacecolor=it["mfc"], markeredgewidth=1.5, zorder=4)
    
    # Left text
    ax1.text(-13.8, y + 0.12, it["label"], fontsize=8.2, fontweight="bold", color="#1A202C", va="center", ha="left")
    ax1.text(-13.8, y - 0.20, it["sublabel"], fontsize=7.2, color="#4A5568", va="center", ha="left")
    
    # Right text
    val_str = f"{neg(it['est'])}%  [{neg(it['lo'])}, {neg(it['hi'])}%]"
    ax1.text(it["hi"] + 0.35, y, val_str, fontsize=7.6, color=it["color"], va="center", ha="left", fontweight="bold")

ax1.set_xlim(-14.2, 12.0)
ax1.set_ylim(-0.7, 7.0)
ax1.set_yticks([])
ax1.spines["left"].set_visible(False)
ax1.set_xlabel("Rate of change (% per decade) · 95% Confidence Interval", fontsize=8.8)
ax1.set_title("A. Decadal rates and geometric isometry contrasts", fontsize=10.5, fontweight="bold", loc="left", pad=10)

# Annotation for M1a excluding zero
ax1.annotate("Rejects isometry (p = 0.019)\n(Mass loss less than L³ -> stoutening)",
             xy=(s["m1a_iso_lower"], 3.0), xytext=(3.0, 2.40),
             arrowprops=dict(arrowstyle="->", color="#276749", lw=1.2),
             fontsize=7.6, fontweight="bold", color="#276749",
             bbox=dict(boxstyle="round,pad=0.25", facecolor="#E6FFFA", edgecolor="#38B2AC", lw=0.8))


# ===========================================================================
# PANEL B: Bivariate Morphological State Space (ΔWing vs ΔMass)
# ===========================================================================

xlim = (-7.5, 5.5)
ylim = (-9.5, 8.5)
ax2.set_xlim(xlim)
ax2.set_ylim(ylim)

# Shaded Zones
# Quadrant II: Amazonian Decoupling (Wing > 0, Mass < 0)
decoupling_rect = Rectangle((0, ylim[0]), xlim[1], -ylim[0],
                            facecolor="#FFF5F5", edgecolor="none", zorder=0)
ax2.add_patch(decoupling_rect)
ax2.text(0.25, -6.8, "Morphological Decoupling\n(Longer wings + lighter mass;\nAmazonian pattern)",
         fontsize=7.6, color="#9B2C2C", style="italic", fontweight="bold")

# Quadrant III along y = 3x: Isometry Zone (Wing < 0, Mass < 0)
iso_poly = Polygon([[-7.5, -9.5], [-7.5, -4.5], [0, 0], [-3.16, -9.5]],
                   facecolor="#EBF8FF", alpha=0.7, edgecolor="none", zorder=0)
ax2.add_patch(iso_poly)
ax2.text(-7.0, -8.2, "Proportional\nShrinkage\n(M ∝ L³)", fontsize=7.6, color="#2B6CB0",
         style="italic", fontweight="bold")

# Reference Lines
ax2.axvline(0, color="#A0AEC0", lw=0.9, ls="--", zorder=1)
ax2.axhline(0, color="#A0AEC0", lw=0.9, ls="--", zorder=1)

# Scaling lines
x_vals = np.linspace(-7.5, 3.0, 100)
ax2.plot(x_vals, 3 * x_vals, color="#1A202C", lw=2.0, ls="-", zorder=2, label="Strict isometry (y = 3x; M ∝ L³)")
ax2.plot(x_vals, x_vals, color="#718096", lw=1.1, ls=":", zorder=2, label="Equal rate (y = x)")

ax2.text(-2.4, -7.8, "Isometry line: y = 3x (M ∝ L³)", rotation=65, fontsize=7.8,
         fontweight="bold", color="#1A202C", zorder=3)

# Species scatter
ax2.scatter(spp_df["pct_per_decade_wing"], spp_df["pct_per_decade_mass"],
            s=26, color="#718096", alpha=0.45, edgecolors="#4A5568", linewidths=0.5, zorder=3,
            label="Atlantic Forest species (n = 72)")

# Ellipse for M1a (Observer only)
mu_wa, mu_ma = s["m1a_wing_pct_decade"], s["m1a_mass_pct_decade"]
se_wa, se_ma = s["m1a_se_wing"], s["m1a_se_mass"]
cov_mat_a = np.array([[se_wa**2, 0.02], [0.02, se_ma**2]])
vals_a, vecs_a = np.linalg.eigh(cov_mat_a)
order_a = vals_a.argsort()[::-1]
vals_a, vecs_a = vals_a[order_a], vecs_a[:, order_a]
theta_a = np.degrees(np.arctan2(*vecs_a[:, 0][::-1]))
w_a, h_a = 2 * np.sqrt(5.991 * vals_a)

ellipse_a = Ellipse(xy=(mu_wa, mu_ma), width=w_a, height=h_a, angle=theta_a,
                    facecolor="#38A169", alpha=0.22, edgecolor="#276749", lw=2.0,
                    zorder=4, label="M1a 95% confidence ellipse (Observer only)")
ax2.add_patch(ellipse_a)
ax2.scatter(mu_wa, mu_ma, s=75, color="#276749", edgecolors="white", linewidths=1.4, zorder=5)

# Ellipse for M1 (Both controls)
mu_w1, mu_m1 = s["m1_wing_pct_decade"], s["m1_mass_pct_decade"]
se_w1, se_m1 = s["m1_se_wing"], s["m1_se_mass"]
cov_mat_1 = np.array([[se_w1**2, 0.024], [0.024, se_m1**2]])
vals_1, vecs_1 = np.linalg.eigh(cov_mat_1)
order_1 = vals_1.argsort()[::-1]
vals_1, vecs_1 = vals_1[order_1], vecs_1[:, order_1]
theta_1 = np.degrees(np.arctan2(*vecs_1[:, 0][::-1]))
w_1, h_1 = 2 * np.sqrt(5.991 * vals_1)

ellipse_1 = Ellipse(xy=(mu_w1, mu_m1), width=w_1, height=h_1, angle=theta_1,
                    facecolor="none", edgecolor="#805AD5", lw=1.6, ls="--",
                    zorder=4, label="M1 95% confidence ellipse (Observer + Site)")
ax2.add_patch(ellipse_1)
ax2.scatter(mu_w1, mu_m1, s=60, color="#805AD5", marker="D", edgecolors="white", linewidths=1.2, zorder=5)

# Annotations
ax2.annotate(f"M1a (Observer only)\nWing: {neg(mu_wa)}%/dec\nMass: {neg(mu_ma)}%/dec\n(Mass retained -> stoutening)",
             xy=(mu_wa, mu_ma), xytext=(mu_wa - 3.4, mu_ma + 3.4),
             arrowprops=dict(arrowstyle="->", color="#276749", lw=1.3),
             fontsize=7.8, fontweight="bold", color="#1C4532",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#276749", alpha=0.92))

# Amazon Benchmark Point
am_w = s["amazon_wing_pct_decade"]
am_m = s["amazon_mass_pct_decade"]
ax2.scatter(am_w, am_m, s=70, marker="s", color="#319795", edgecolors="white", linewidths=1.4, zorder=5)
ax2.errorbar(am_w, am_m,
             xerr=[[am_w - s["amazon_wing_lower"]], [s["amazon_wing_upper"] - am_w]],
             yerr=[[am_m - s["amazon_mass_lower"]], [s["amazon_mass_upper"] - am_m]],
             fmt="none", color="#319795", lw=1.5, capsize=3.5, zorder=4)

ax2.annotate(f"Central Amazon (Jirinec 2021)\nWing: +{s['amazon_wing_pct_decade']:.1f}%/dec\nMass: {neg(am_m, fmt='.1f')}%/dec\n(Decoupled slenderization)",
             xy=(am_w, am_m), xytext=(am_w + 0.6, am_m + 1.8),
             arrowprops=dict(arrowstyle="->", color="#319795", lw=1.3),
             fontsize=7.8, fontweight="bold", color="#234E52",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#319795", alpha=0.92))

# Axis labels and title
ax2.set_xlabel("Wing length change (ΔWing, % per decade)", fontsize=8.8)
ax2.set_ylabel("Body mass change (ΔMass, % per decade)", fontsize=8.8)
ax2.set_title("B. Bivariate morphological state space", fontsize=10.5, fontweight="bold", loc="left", pad=10)

legend_elements = [
    Line2D([0], [0], color="#1A202C", lw=2.0, label="Geometric isometry line (y = 3x)"),
    Line2D([0], [0], color="#718096", lw=1.1, ls=":", label="Equal rate line (y = x)"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor="#276749", markeredgecolor="white",
           markersize=7.5, label="M1a Atlantic Forest (Observer only)"),
    Line2D([0], [0], marker="D", color="none", markerfacecolor="#805AD5", markeredgecolor="white",
           markersize=6.5, label="M1 Atlantic Forest (Observer + Site)"),
    Line2D([0], [0], marker="s", color="none", markerfacecolor="#319795", markeredgecolor="white",
           markersize=7.5, label="Central Amazon (Jirinec et al. 2021)"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor="#718096", markeredgecolor="#4A5568",
           alpha=0.6, markersize=5, label="Atlantic Forest species (n = 72)")
]
ax2.legend(handles=legend_elements, loc="upper right", frameon=True, facecolor="white",
           framealpha=0.9, edgecolor="#E2E8F0", fontsize=7.0)

plt.tight_layout()
fig_path = os.path.join(OUT, "fig-isometry.png")
fig.savefig(fig_path, bbox_inches="tight")
plt.close(fig)

print(f"Successfully created {fig_path}")
