#!/usr/bin/env python3
"""Regenerate manuscript figures from the source data.

Figures produced here (saved to Manuscript/images/):
  fig-map.png         Fig 1  sampling localities over a South America basemap
  fig-wingtrend.png   Fig 2  wing length vs year with the model-estimated trend
  fig-trends.png      Fig 3  per-species wing & bill trajectories (early vs late)
  fig-variability.png Fig 4  lnCVR forest plots, two panels (wing & bill)

The arthropod figure (fig-arthropods.png) is exported from the Rmd
(`export_fig_arthropods` chunk) because it is the fitted brms model.

Basemap geometry is read from Manuscript/data/south_america.geojson (committed,
offline). Model estimates below come from atlantic_parallel.R (rubin_summary).
"""
import json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.join(HERE, "images")
DATA = os.path.join(HERE, "data", "south_america.geojson")
# Canonical analytical sample (89 species / 15,332 records), exported from
# passer90.rda by Analysis/update_descriptive_stats.R so the figures use the
# exact same build as the models and the quoted text.
CSV  = os.path.join(HERE, "..", "Analysis", "data", "derived", "passer90_export.csv")

# --- model-estimated wing~year effect (Rubin-pooled; atlantic_parallel.R) ---
WING_YR_BETA, WING_YR_LO, WING_YR_HI = -0.7294, -1.0854, -0.3733   # mm per SD-year
# --- canonical lnCVR meta-analytic means (metafor::rma REML; descriptive_summary.rds) ---
# Used for the summary line so the figure matches the values quoted in the text.
LNCVR = {"cwl": (0.241, 0.104, 0.378), "Bill_width.mm.": (-0.073, -0.276, 0.131)}

plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 150})
DEC, INC = "#CC6666", "#5B7FB5"

# --- canonical analytical sample (already filtered and name-reconciled) -----
f = pd.read_csv(CSV, low_memory=False)
f["Binomial"] = f["Binomial"].str.replace(" ", "_")
print("species", f.Binomial.nunique(), "records", len(f))

q = f[(f.Year <= 2006) | (f.Year >= 2013)].copy()
q["period"] = np.where(q.Year <= 2006, "1990–2006", "2013–2018")

# ---------------------------------------------------------------------------
# Basemap helper: draw South America polygons from the committed GeoJSON.
# ---------------------------------------------------------------------------
def draw_basemap(ax):
    try:
        gj = json.load(open(DATA))
    except Exception as e:
        print("basemap skipped:", e); return
    def rings(geom):
        if geom["type"] == "Polygon":
            return [geom["coordinates"]]
        if geom["type"] == "MultiPolygon":
            return geom["coordinates"]
        return []
    for feat in gj["features"]:
        for poly in rings(feat["geometry"]):
            ext = np.array(poly[0])
            ax.fill(ext[:, 0], ext[:, 1], facecolor="#ECECEC", edgecolor="#9A9A9A",
                    linewidth=0.4, zorder=0)

# ---------------------------------------------------------------------------
# FIG 1 - sampling localities over the basemap, by period
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(8, 4.4))
gs = gridspec.GridSpec(1, 2, wspace=0.18)
vmin, vmax = q.cwl.quantile(.02), q.cwl.quantile(.98)
for i, p in enumerate(["1990–2006", "2013–2018"]):
    ax = fig.add_subplot(gs[0, i]); d = q[q.period == p].dropna(subset=["cwl"])
    draw_basemap(ax)
    sc = ax.scatter(d.Longitude_decimal_degrees, d.Latitude_decimal_degrees, c=d.cwl,
                    cmap="plasma", s=9, alpha=0.6, edgecolors="none", vmin=vmin, vmax=vmax, zorder=2)
    ax.set_title(p, fontsize=10); ax.set_xlabel("Longitude (°)")
    if i == 0: ax.set_ylabel("Latitude (°)")
    ax.set_xlim(-60, -32); ax.set_ylim(-34, -3); ax.set_aspect("equal", "box")
cb = fig.colorbar(sc, ax=fig.axes, fraction=0.025, pad=0.02); cb.set_label("Wing length (mm)")
fig.suptitle(f"Sampling localities of {f.Binomial.nunique()} Atlantic Forest passerine species", fontsize=11, y=0.99)
fig.savefig(f"{OUT}/fig-map.png", bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------------------
# FIG 2 - wing length vs year, with the model-estimated trend (NEW)
# ---------------------------------------------------------------------------
w = f.dropna(subset=["cwl"])
ybar = w.cwl.mean(); ymu = w.Year.mean(); ysd = w.Year.std()
yr_grid = np.linspace(w.Year.min(), w.Year.max(), 100)
z = (yr_grid - ymu) / ysd
fit = ybar + WING_YR_BETA * z
lo  = ybar + WING_YR_LO * z
hi  = ybar + WING_YR_HI * z
fig, ax = plt.subplots(figsize=(6.4, 4.6))
# transparent points (small x-jitter so within-year records don't stack on a line)
jit = (np.random.default_rng(1).random(len(w)) - 0.5) * 0.7
ax.scatter(w.Year + jit, w.cwl, s=6, alpha=0.08, color="#3A3A3A", edgecolors="none", zorder=1)
ax.fill_between(yr_grid, lo, hi, color="#5B7FB5", alpha=0.35, zorder=2, label="95% CI")
ax.plot(yr_grid, fit, color="#2E5090", lw=2, zorder=3, label="Model-estimated trend")
ax.set_xlabel("Year"); ax.set_ylabel("Wing length (mm)")
ax.set_ylim(w.cwl.quantile(.005), w.cwl.quantile(.995))
ax.set_title("Wing length decline, 1990–2018\n"
             f"year effect β = {WING_YR_BETA:.2f} mm per SD-year "
             f"[{WING_YR_LO:.2f}, {WING_YR_HI:.2f}]".replace("-", "−"), fontsize=10)
ax.legend(loc="upper right", frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(f"{OUT}/fig-wingtrend.png", bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------------------
# FIG 3 - per-species wing & bill trajectories, early vs late
# ---------------------------------------------------------------------------
def per_species(metric):
    g = q.dropna(subset=[metric]).groupby(["Binomial", "period"])[metric].mean()
    return g.unstack().dropna()
fig, axes = plt.subplots(1, 2, figsize=(8, 4.3))
for ax, (metric, lab) in zip(axes, [("cwl", "Wing length"), ("Bill_width.mm.", "Bill width")]):
    wmat = per_species(metric)
    for sp, row in wmat.iterrows():
        e, l = row["1990–2006"], row["2013–2018"]
        ax.plot([0, 1], [e, l], color=(DEC if l < e else INC), alpha=0.6, lw=1)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["1990–2006", "2013–2018"])
    ax.set_ylabel(f"{lab} (mm)"); ax.set_title(lab, fontsize=10)
fig.legend([Line2D([0], [0], color=DEC), Line2D([0], [0], color=INC)],
           ["Decrease", "Increase"], loc="lower center", ncol=2, frameon=False,
           bbox_to_anchor=(0.5, -0.04))
fig.tight_layout(); fig.savefig(f"{OUT}/fig-trends.png", bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------------------
# FIG 4 - lnCVR forest plots, TWO panels (wing & bill); corrected SE
# ---------------------------------------------------------------------------
def lncvr_table(metric):
    g = q.dropna(subset=[metric]).groupby(["Binomial", "period"])[metric].agg(m="mean", s="std", n="count").reset_index()
    piv = g.pivot(index="Binomial", columns="period").dropna()
    me, mc = piv[("m", "2013–2018")], piv[("m", "1990–2006")]
    se, sc = piv[("s", "2013–2018")], piv[("s", "1990–2006")]
    ne, nc = piv[("n", "2013–2018")], piv[("n", "1990–2006")]
    def corr(period):
        mm = np.log(g[g.period == period]["m"]); ss = np.log(g[g.period == period]["s"])
        ok = np.isfinite(mm) & np.isfinite(ss); return np.corrcoef(mm[ok], ss[ok])[0, 1]
    ce, cc = corr("2013–2018"), corr("1990–2006")
    lncvr = np.log((se / me) / (sc / mc)) + 1 / (2 * (ne - 1)) - 1 / (2 * (nc - 1))
    def s2(mc, sc, nc, corc, me, se, ne, core):
        ac = sc**2 / (nc * mc**2); bc = 1 / (2 * (nc - 1)); cc_ = 2 * corc * np.sqrt(ac * bc)
        ae = se**2 / (ne * me**2); be = 1 / (2 * (ne - 1)); ce_ = 2 * core * np.sqrt(ae * be)
        return ac + bc - cc_ + ae + be - ce_
    var = s2(mc, sc, nc, cc, me, se, ne, ce)
    t = pd.DataFrame({"sp": piv.index, "lncvr": lncvr.values, "var": var.values}).dropna()
    t["se"] = np.sqrt(t["var"]); return t.sort_values("lncvr")

def re_mean(t):
    wts = 1 / t["var"]; mu = (wts * t["lncvr"]).sum() / wts.sum()
    Q = (wts * (t["lncvr"] - mu) ** 2).sum(); k = len(t)
    C = wts.sum() - (wts ** 2).sum() / wts.sum()
    tau2 = max(0, (Q - (k - 1)) / C); wr = 1 / (t["var"] + tau2)
    mur = (wr * t["lncvr"]).sum() / wr.sum(); se = np.sqrt(1 / wr.sum())
    return mur, mur - 1.96 * se, mur + 1.96 * se

fig, axes = plt.subplots(1, 2, figsize=(9, 8))
for ax, (metric, lab) in zip(axes, [("cwl", "Wing length"), ("Bill_width.mm.", "Bill width")]):
    t = lncvr_table(metric); y = np.arange(len(t))
    mu, mlo, mhi = LNCVR[metric]                       # canonical metafor values (match the text)
    chk = re_mean(t)                                   # Python re-computation, for a sanity check
    ax.errorbar(t["lncvr"], y, xerr=1.96 * t["se"], fmt="o", ms=2.5, lw=0.6,
                color="#333333", ecolor="#AAAAAA", capsize=1.2)
    ax.axvline(0, ls="--", color="grey", lw=0.8)
    ax.axvspan(mlo, mhi, color="#5B7FB5", alpha=0.18); ax.axvline(mu, color="#2E5090", lw=1.3)
    ax.set_yticks(y); ax.set_yticklabels([s.replace("_", " ") for s in t["sp"]], fontsize=4.5, style="italic")
    ax.set_xlabel("lnCVR (late vs early) · 95% CI")
    ax.set_title(f"{lab}\nmeta-analytic mean = {mu:.2f} [{mlo:.2f}, {mhi:.2f}]", fontsize=9)
    print(f"{lab}: metafor mean {mu} [{mlo}, {mhi}]; python re-check {chk[0]:.3f} [{chk[1]:.3f}, {chk[2]:.3f}] (k={len(t)})")
fig.tight_layout(); fig.savefig(f"{OUT}/fig-variability.png", bbox_inches="tight"); plt.close(fig)
print("figures written to", OUT)
