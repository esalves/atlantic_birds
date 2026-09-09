#!/usr/bin/env python3
"""Regenerate manuscript figures from the saved analysis outputs.

Figures produced here (saved to Manuscript/images/):
  fig-map.png                 Fig 1  sampling sites, early vs late period; sites sampled in
                                     BOTH periods drawn distinctly from period-unique sites,
                                     N records / sites / species per panel (revision 2026-09)
  fig-wingtrend.png           Fig 2  within-species(x sex)-centred wing length vs year, annual
                                     means, the controlled M3 fit (glmmTMB phylogenetic mixed
                                     model, Rubin-pooled over trees, when available; else lme4
                                     fast tier) against the uncontrolled M0 baseline, plus a
                                     lighter lme4 Tier-1 M3 comparison line, inset of records/yr
  fig-species-slopes.png      new    species-specific year slopes (caterpillar) from the
                                     controlled M3 model (glmmTMB phylogenetic tier when
                                     available), with the M0 and lme4 Tier-1 M3 slopes shown
                                     as comparison layers
  fig-trends.png              Fig 3  per-species wing & bill trajectories (early vs late) [unchanged]
  fig-variability.png         Fig 4  lnCVR forest plots, two panels (wing & bill)         [unchanged]
  fig-s-records-by-source.png Supp.  wing records per year stacked by contributor, with the
                                     number of contributors per year
  fig-s-wingcol.png           Supp.  which wing column supplied the coalesced wing length,
                                     by year (measurement-protocol proxy)

Inputs
------
All inputs are plain CSV files written by Analysis/scripts/make_figures.R (mode
"export") into Analysis/output/figure_data/, flattened from the R objects that
carry the numbers quoted in the manuscript:
  fig_records.csv                  <- data/derived/passer90.rda (live-only, known-sex, 73 spp)
  fig_map_sites.csv, fig_scalars   <- output/audit_sites.rds, audit_sources.rds, effect_scale.rds (P0)
  fig_contributor*.csv, fig_records_per_contributor_year.csv, fig_wingcol_by_year.csv (P0)
  fig_controlled_before_after.csv        <- output/controlled_wing_results.rds (P1, lme4 fast tier)
  fig_species_slopes.csv                 <- output/controlled_wing_species_slopes.rds (P1)
  fig_controlled_before_after_phylo.csv  <- output/controlled_wing_phylo_results.rds (P1, glmmTMB
  fig_species_slopes_phylo.csv           <- output/controlled_wing_phylo_species_slopes.rds  phylo tier, optional)
This script first tries to refresh those CSVs by calling
`Rscript Analysis/scripts/make_figures.R export` (pass --no-export to skip); if
Rscript is unavailable it uses the CSVs already on disk. Figures whose inputs
are missing are skipped with a warning, never drawn from placeholder numbers.

Figs 3 and 4 are unchanged from the submitted version and still read
Analysis/data/derived/passer90_export.csv (written by update_descriptive_stats.R).
Basemap geometry is read from Manuscript/data/south_america.geojson (committed, offline).

Tier priority for the M3 fit in Fig 2 and the slopes in fig-species-slopes
(2026-09 revision, GLMMTMB_ENGINE.md): lme4 REML fast tier (Tier 1 of
REVISION_PLAN.md Phase 1; Wald intervals, no phylogenetic term) is the base
layer and is always drawn as a lighter comparison line/markers once a higher
tier is available. When controlled_wing_phylo_results.rds /
controlled_wing_phylo_species_slopes.rds are present (glmmTMB propto,
Rubin-pooled over the same 50 trees as the published brms model — see
Williams et al. 2025 and glmmtmb_validation_wing.R), those become the primary
fitted line / markers, labelled "glmmTMB phylogenetic mixed model, N trees,
Rubin-pooled" with N read from the file (never hard-coded, so a partial dev
run is never mislabelled as the full run). When controlled_wing_brms_results.rds
(the Tier 2-3 Bayesian cross-check) is also present it takes priority over both
and the panels are relabelled again accordingly.

Run:  python3 make_figures.py            (from Manuscript/ or the repo root)
      python3 make_figures.py --no-export
Environment: python 3.11, matplotlib 3.11, pandas 3.0, numpy 2.x (requirements.txt).
"""
import json, os, sys, shutil, subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.join(HERE, "images")
DATA = os.path.join(HERE, "data", "south_america.geojson")
ANALYSIS = os.path.normpath(os.path.join(HERE, "..", "Analysis"))
FD   = os.path.join(ANALYSIS, "output", "figure_data")          # CSVs from make_figures.R export
# Canonical analytical sample (live-only: 73 species / 12,571 records), exported from
# passer90.rda by Analysis/scripts/update_descriptive_stats.R; still used by Figs 3-4.
CSV  = os.path.join(ANALYSIS, "data", "derived", "passer90_export.csv")

# --- canonical lnCVR meta-analytic means (metafor::rma REML; descriptive_summary.rds) ---
# Used for the Fig 4 summary line so the figure matches the values quoted in the text.
LNCVR = {"cwl": (0.183, 0.033, 0.334), "Bill_width.mm.": (0.107, -0.196, 0.411)}

plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 150})
DEC, INC = "#CC6666", "#5B7FB5"          # decrease / increase (house colours)
SHARED, UNIQUE = "#B2182B", "#4393C3"    # Fig 1: site sampled in both periods / period-unique
EARLY_MAX, LATE_MIN = 2006, 2013         # quartile periods used throughout the manuscript
PERIODS = [("early", f"1995–{EARLY_MAX}"), ("late", f"{LATE_MIN}–2018")]

def neg(s):
    """Typographic minus for figure text."""
    return str(s).replace("-", "−")

# ---------------------------------------------------------------------------
# 0. Refresh the CSV inputs from the R objects (Rscript make_figures.R export)
# ---------------------------------------------------------------------------
if "--no-export" not in sys.argv:
    rscript = shutil.which("Rscript")
    if rscript:
        r = subprocess.run([rscript, os.path.join(ANALYSIS, "scripts", "make_figures.R"), "export"],
                           cwd=HERE, capture_output=True, text=True)
        print(r.stderr.strip() or r.stdout.strip())
        if r.returncode != 0:
            print("WARNING: R export step failed; using the CSVs already in", FD)
    else:
        print("Rscript not found; using the CSVs already in", FD)

def fd_csv(name, required=True, quiet=False):
    p = os.path.join(FD, name)
    if not os.path.exists(p):
        if required:
            sys.exit(f"missing {p}: run `Rscript Analysis/scripts/make_figures.R export` first")
        if not quiet:
            print(f"WARNING: {name} not found - figures that need it are skipped")
        return None
    return pd.read_csv(p, low_memory=False)

rec  = fd_csv("fig_records.csv")
scal = fd_csv("fig_scalars.csv").set_index("key")["value"]
def S(key, cast=float, default=None):
    return cast(scal[key]) if key in scal.index else default
YEAR_CENTRE = S("year_centre_scaling")      # 2009.487: scaled_yr = (Year - centre) / sd
YEAR_SD     = S("year_sd_scaling")          # 5.0199  (the SD the coefficients are per)
N_SPP       = rec.Binomial.nunique()
print(f"records {len(rec)}  species {N_SPP}  wing records {rec.cwl.notna().sum()}  "
      f"scaled_yr centre {YEAR_CENTRE:.3f} sd {YEAR_SD:.4f}")

# numbers printed on the figures, written back for the manuscript to quote
panel_numbers = []
def note(figure, item, value):
    panel_numbers.append({"figure": figure, "item": item, "value": value})

# --- canonical sample for Figs 3-4 (unchanged code path) ----------------------------
f = pd.read_csv(CSV, low_memory=False)
f["Binomial"] = f["Binomial"].str.replace(" ", "_")
if len(f) != len(rec) or f.Binomial.nunique() != N_SPP:
    print(f"WARNING: passer90_export.csv ({len(f)} rows, {f.Binomial.nunique()} spp) differs from "
          f"fig_records.csv ({len(rec)} rows, {N_SPP} spp); re-run update_descriptive_stats.R")
q = f[(f.Year <= EARLY_MAX) | (f.Year >= LATE_MIN)].copy()
q["period"] = np.where(q.Year <= EARLY_MAX, PERIODS[0][1], PERIODS[1][1])

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
# FIG 1 - sampling sites by period: shared vs period-unique, N per panel
# ---------------------------------------------------------------------------
# One row per coordinate site (exact lon/lat pair) with wing records in the quartile
# periods (audit_sites.rds$map_sites, P0). `shared` = records in both periods.
sites = fd_csv("fig_map_sites.csv").dropna(subset=["Longitude_decimal_degrees", "Latitude_decimal_degrees"])
wq = rec[rec.cwl.notna() & ((rec.Year <= EARLY_MAX) | (rec.Year >= LATE_MIN))].copy()
wq["period"] = np.where(wq.Year <= EARLY_MAX, "early", "late")
n_shared_sites = int(sites["shared"].sum())
fig = plt.figure(figsize=(8.4, 4.9))
gs = gridspec.GridSpec(1, 2, wspace=0.12)
smax = sites[["n_early", "n_late"]].to_numpy().max()
def msize(n):  # area proportional to records
    return 6 + 140 * np.sqrt(n / smax)
for i, (per, lab) in enumerate(PERIODS):
    ax = fig.add_subplot(gs[0, i])
    draw_basemap(ax)
    d = sites[sites[f"n_{per}"] > 0]
    u, s = d[~d["shared"]], d[d["shared"]]
    ax.scatter(u.Longitude_decimal_degrees, u.Latitude_decimal_degrees, s=msize(u[f"n_{per}"]),
               facecolors=UNIQUE, edgecolors="white", linewidths=0.4, alpha=0.75, zorder=2)
    ax.scatter(s.Longitude_decimal_degrees, s.Latitude_decimal_degrees, s=msize(s[f"n_{per}"]) * 1.15,
               marker="^", facecolors=SHARED, edgecolors="black", linewidths=0.7, alpha=0.95, zorder=3)
    n_rec = int(d[f"n_{per}"].sum()); n_site = len(d); n_spp = wq.loc[wq.period == per, "Binomial"].nunique()
    n_contrib = wq.loc[wq.period == per, "Main_researcher"].nunique()
    for k, v in [("records", n_rec), ("coordinate sites", n_site), ("shared sites", int(d["shared"].sum())),
                 ("species", n_spp), ("contributors", n_contrib)]:
        note("fig-map", f"{per}: {k}", v)
    ax.text(0.03, 0.97, f"N = {n_rec:,} wing records\n{n_site} sites ({int(d['shared'].sum())} sampled in both periods)\n"
            f"{n_spp} species, {n_contrib} contributors", transform=ax.transAxes, fontsize=7.5,
            va="top", ha="left", zorder=5,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="none", alpha=0.85))
    ax.set_title(f"{lab}", fontsize=10); ax.set_xlabel("Longitude (°)")
    if i == 0: ax.set_ylabel("Latitude (°)")
    ax.set_xlim(-60, -32); ax.set_ylim(-34, -3); ax.set_aspect("equal", "box")
handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=UNIQUE, markeredgecolor="white", markersize=7,
                  label="site sampled in one period only"),
           Line2D([0], [0], marker="^", color="none", markerfacecolor=SHARED, markeredgecolor="black", markersize=8,
                  label=f"site sampled in both periods (n = {n_shared_sites})"),
           Line2D([0], [0], marker="o", color="none", markerfacecolor="#999999", markersize=4, label="marker area ∝ records")]
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.02))
mun_sh, mun_n = S("municipality_wing.n_shared", int), S("municipality_wing.n_units", int)
loc_sh, loc_n = S("locality_wing.n_shared", int), S("locality_wing.n_units", int)
fig.suptitle(f"Sampling sites of {N_SPP} Atlantic Forest passerine species, early vs late period\n"
             f"coordinate sites sampled in both periods: {n_shared_sites} of {len(sites)}; named localities {loc_sh} of {loc_n}; "
             f"municipalities {mun_sh} of {mun_n}", fontsize=9.5, y=1.0)
fig.savefig(f"{OUT}/fig-map.png", bbox_inches="tight"); plt.close(fig)
print(f"fig-map: sites early {int((sites.n_early > 0).sum())}, late {int((sites.n_late > 0).sum())}, shared {n_shared_sites}")

# ---------------------------------------------------------------------------
# FIG 2 - within-species-centred wing length vs year, with the controlled (M3) fit
# ---------------------------------------------------------------------------
# Each record is expressed as its deviation from the species x sex mean wing length
# (Sex is a fixed effect in every model), so species with different absolute sizes
# are put on a common scale and the pooled-species artefact the referee objected to
# (§2.3) is removed. The fitted lines are the year coefficients of before_after
# tables built the same way at every tier: M0 = the published structure (species +
# phylogeny intercepts/slopes, Sex, latitude) and M3 = M0 + contributor,
# municipality, wing-column proxy, longitude, altitude, season, moult and individual
# terms. A line is beta * (Year - centre) / SD(year), centred on the sample mean year,
# i.e. the marginal year effect holding every other term fixed.
#
# Tier priority (2026-09 revision, GLMMTMB_ENGINE.md): lme4 Tier 1 is always read
# and, once a higher tier is available, kept only as a lighter non-phylogenetic
# comparison line; the glmmTMB phylogenetic tier (propto, Rubin-pooled over trees;
# Williams et al. 2025, validated by glmmtmb_validation_wing.R) is promoted to the
# primary fit when present; the brms Tier 2-3 Bayesian cross-check, if copied back
# from Totoro, is promoted above that in turn.
ba_lme4  = fd_csv("fig_controlled_before_after.csv", required=False)
ba_phylo = fd_csv("fig_controlled_before_after_phylo.csv", required=False, quiet=True) if ba_lme4 is not None else None
ba_brms  = fd_csv("fig_controlled_before_after_brms.csv", required=False, quiet=True) if ba_lme4 is not None else None
ba = ba_lme4
ba_compare = None   # lighter comparison line drawn alongside the primary fit, when available
tier_label = "lme4 fast tier (Tier 1, REML, Wald 95% CI; phylogenetic fit pending)"
if ba_phylo is not None:
    ba_compare = ba_lme4
    n_trees_phylo = S("phylo_n_trees", int, default=None)
    tier_label = (f"glmmTMB phylogenetic mixed model, {n_trees_phylo} trees, Rubin-pooled (Wald 95% CI)"
                  if n_trees_phylo is not None else "glmmTMB phylogenetic mixed model, Rubin-pooled (Wald 95% CI)")
    ba = ba_phylo
if ba_brms is not None:
    ba, tier_label = ba_brms, "brms, Rubin-pooled across phylogenetic trees (95% CrI)"
def year_term(df, model, term="scaled_yr"):
    if df is None:
        return None
    r = df[(df.model == model) & (df.term == term)]
    if len(r) != 1:
        return None
    r = r.iloc[0]
    return dict(est=float(r.estimate), lo=float(r.ci_lo), hi=float(r.ci_hi), N=int(r.N),
                mmd=float(r.mm_per_decade), mmd_lo=float(r.mm_per_decade_lo), mmd_hi=float(r.mm_per_decade_hi))

w = rec.dropna(subset=["cwl"]).copy()
w["dev"] = w.cwl - w.groupby(["Binomial", "Sex"]).cwl.transform("mean")
ann = w.groupby("Year").dev.agg(mean="mean", sd="std", n="count").reset_index()
ann["se"] = ann["sd"] / np.sqrt(ann["n"])
yr_grid = np.linspace(w.Year.min(), w.Year.max(), 100)
z = (yr_grid - YEAR_CENTRE) / YEAR_SD

fig = plt.figure(figsize=(6.8, 7.4))
gs = gridspec.GridSpec(2, 1, height_ratios=[1.35, 1], hspace=0.08)
ax = fig.add_subplot(gs[0]); axz = fig.add_subplot(gs[1], sharex=ax)
jit = (np.random.default_rng(1).random(len(w)) - 0.5) * 0.7
ax.scatter(w.Year + jit, w.dev, s=5, alpha=0.07, color="#3A3A3A", edgecolors="none", zorder=1)
m0 = m3 = m3_compare = None
if ba is not None:
    m0, m3 = year_term(ba, "M0"), year_term(ba, "M3")
if ba_compare is not None:
    m3_compare = year_term(ba_compare, "M3")
for a in (ax, axz):
    a.errorbar(ann.Year, ann["mean"], yerr=1.96 * ann.se, fmt="o", ms=3.2, color="black", ecolor="#555555",
               elinewidth=0.8, capsize=1.5, zorder=4, label="annual mean ± 95% CI")
    a.axhline(0, color="#BBBBBB", lw=0.6, zorder=0)
    if m0:
        a.plot(yr_grid, m0["est"] * z, color="#777777", lw=1.6, ls="--", zorder=3,
               label=f"M0 baseline: {neg(f'{m0['mmd']:.2f}')} mm/decade [{neg(f'{m0['mmd_lo']:.2f}')}, {neg(f'{m0['mmd_hi']:.2f}')}]")
    if m3_compare:
        a.plot(yr_grid, m3_compare["est"] * z, color="#2E5090", lw=1.0, ls=":", alpha=0.55, zorder=2,
               label=f"lme4 Tier-1 M3 (comparison): {neg(f'{m3_compare['mmd']:.2f}')} mm/decade "
                     f"[{neg(f'{m3_compare['mmd_lo']:.2f}')}, {neg(f'{m3_compare['mmd_hi']:.2f}')}]")
    if m3:
        a.fill_between(yr_grid, m3["lo"] * z, m3["hi"] * z, color=INC, alpha=0.30, zorder=2, linewidth=0)
        a.plot(yr_grid, m3["est"] * z, color="#2E5090", lw=2.2, zorder=3,
               label=f"M3 full controls: {neg(f'{m3['mmd']:.2f}')} mm/decade [{neg(f'{m3['mmd_lo']:.2f}')}, {neg(f'{m3['mmd_hi']:.2f}')}]")
if m3:
    for k in ("est", "lo", "hi", "mmd", "mmd_lo", "mmd_hi", "N"):
        note("fig-wingtrend", f"M3 {k}", m3[k]); note("fig-wingtrend", f"M0 {k}", m0[k])
    note("fig-wingtrend", "M3 tier", tier_label)
    if m3_compare:
        for k in ("est", "lo", "hi", "mmd", "mmd_lo", "mmd_hi", "N"):
            note("fig-wingtrend", f"M3_lme4_compare {k}", m3_compare[k])
else:
    print("WARNING: fig-wingtrend drawn WITHOUT model lines (controlled results missing)")
ax.set_ylim(w.dev.quantile(.005), w.dev.quantile(.995)); ax.tick_params(labelbottom=False)
ax.set_ylabel("Deviation from species × sex\nmean wing length (mm)")
ax.set_title(f"Within-species-centred wing length, 1995–2018 (N = {len(w):,} records, {w.Binomial.nunique()} species)\n"
             f"fitted year effects: {tier_label}", fontsize=9)
ax.text(0.01, 0.98, "a", transform=ax.transAxes, fontweight="bold", va="top")
# inset: wing records per year (top-left, where the early years have few extreme records)
ins = ax.inset_axes([0.06, 0.70, 0.28, 0.27])
ins.bar(ann.Year, ann.n, width=0.8, color="#9E9E9E", edgecolor="none")
ins.set_title("wing records per year", fontsize=6.5, pad=2)
ins.tick_params(labelsize=5.5, length=2, pad=1)
ins.set_xticks([1995, 2000, 2005, 2010, 2015]); ins.spines[["top", "right"]].set_visible(False)
ins.patch.set_alpha(0.85)
# zoom panel: annual means and the fitted lines on the scale of the effect
zl = max(3.0, float(np.nanmax(np.abs(ann["mean"] + np.sign(ann["mean"]) * 1.96 * ann.se))) * 1.05)
axz.set_ylim(-zl, zl); axz.set_xlabel("Year")
axz.set_ylabel("Deviation from species × sex\nmean wing length (mm), zoom")
axz.text(0.01, 0.98, "b", transform=axz.transAxes, fontweight="bold", va="top")
axz.legend(loc="lower left", frameon=False, fontsize=7.2)
fig.savefig(f"{OUT}/fig-wingtrend.png", bbox_inches="tight"); plt.close(fig)
ann.assign(figure="fig-wingtrend").to_csv(os.path.join(FD, "fig_wingtrend_annual_means.csv"), index=False)
if m3:
    print(f"fig-wingtrend: M0 {m0['est']:.3f} [{m0['lo']:.3f}, {m0['hi']:.3f}]  M3 {m3['est']:.3f} [{m3['lo']:.3f}, {m3['hi']:.3f}] "
          f"= {m3['mmd']:.2f} mm/decade, N = {m3['N']}")

# ---------------------------------------------------------------------------
# NEW FIG - species-specific year slopes (caterpillar) from the controlled model
# ---------------------------------------------------------------------------
# slope = fixed scaled_yr + species random slope; interval = +/- 1.96 * sqrt(SE_fixed^2 +
# SD_BLUP^2) — quadrature SE, same convention at every tier (see the source .rds
# $note). Filled markers: interval excludes zero. The M0 (baseline, uncontrolled)
# slope of each species is drawn as a small hollow grey marker, and (once a higher
# tier is primary) the lme4 Tier-1 M3 slope as a light "x" marker, so the shift
# across controls and across tiers is both visible.
#
# Tier priority mirrors Fig 2 (GLMMTMB_ENGINE.md): lme4 Tier 1 -> glmmTMB
# phylogenetic tier (if present) -> brms Tier 2-3 (if present), each promoted over
# the last as the primary caterpillar, with the lower tiers kept as comparison layers.
sl_lme4  = fd_csv("fig_species_slopes.csv", required=False)
sl_phylo = fd_csv("fig_species_slopes_phylo.csv", required=False, quiet=True) if sl_lme4 is not None else None
sl_brms  = fd_csv("fig_species_slopes_brms_M3.csv", required=False, quiet=True) if sl_lme4 is not None else None
sl = sl_lme4
if sl is not None:
    src_label = "lme4 fast tier (Tier 1; phylogenetic fit pending)"
    d3 = sl[sl.model == "M3"].copy()
    d3_lme4_compare = None
    if sl_phylo is not None and {"spp", "mm_per_decade", "mm_per_decade_lo", "mm_per_decade_hi"} <= set(sl_phylo.columns):
        d3_lme4_compare = d3.copy()
        n_trees_phylo = S("phylo_n_trees", int, default=None)
        src_label = (f"glmmTMB phylogenetic mixed model, {n_trees_phylo} trees, Rubin-pooled"
                     if n_trees_phylo is not None else "glmmTMB phylogenetic mixed model, Rubin-pooled")
        d3 = sl_phylo[sl_phylo.model == "M3"].copy()
    if sl_brms is not None and {"spp", "mm_per_decade", "mm_per_decade_lo", "mm_per_decade_hi"} <= set(sl_brms.columns):
        d3, src_label = sl_brms.copy(), "brms, Rubin-pooled across phylogenetic trees"
    d0 = sl[sl.model == "M0"].set_index("spp")
    d3_cmp = d3_lme4_compare.set_index("spp") if d3_lme4_compare is not None else None
    d3 = d3.sort_values("mm_per_decade").reset_index(drop=True)
    y = np.arange(len(d3))
    excl = (d3.mm_per_decade_hi < 0) | (d3.mm_per_decade_lo > 0)
    cols = np.where(d3.mm_per_decade < 0, DEC, INC)
    fig, ax = plt.subplots(figsize=(6.6, 10.5))
    ax.axvline(0, color="grey", lw=0.8, ls="--", zorder=1)
    if m3:
        ax.axvspan(m3["mmd_lo"], m3["mmd_hi"], color=INC, alpha=0.15, zorder=0)
        ax.axvline(m3["mmd"], color="#2E5090", lw=1.3, zorder=1)
    if len(d0):
        ax.scatter(d0.reindex(d3.spp).mm_per_decade, y, s=12, facecolors="none", edgecolors="#8A8A8A",
                   linewidths=0.7, zorder=2, label="M0 baseline (uncontrolled) slope")
    if d3_cmp is not None and len(d3_cmp):
        ax.scatter(d3_cmp.reindex(d3.spp).mm_per_decade, y, s=14, marker="x", color="#2E5090", alpha=0.45,
                   linewidths=0.8, zorder=2, label="lme4 Tier-1 M3 (comparison) slope")
    ax.hlines(y, d3.mm_per_decade_lo, d3.mm_per_decade_hi, colors=cols, lw=0.9, alpha=0.8, zorder=3)
    ax.scatter(d3.mm_per_decade[excl], y[excl], s=18, c=cols[excl], edgecolors=cols[excl], zorder=4)
    ax.scatter(d3.mm_per_decade[~excl], y[~excl], s=18, facecolors="white", edgecolors=cols[~excl], linewidths=1.0, zorder=4)
    ax.set_yticks(y); ax.set_yticklabels([s.replace("_", " ") for s in d3.spp], fontsize=5.6, style="italic")
    ax.set_ylim(-1, len(d3)); ax.set_xlabel("Species year slope, M3 full controls (mm per decade) · 95% interval")
    # x-range covers every point estimate; intervals wider than the axis are clipped and flagged
    XLIM = 10.0
    xlo = min(-XLIM, np.floor(d3.mm_per_decade.min()) - 1); xhi = max(XLIM, np.ceil(d3.mm_per_decade.max()) + 1)
    ax.set_xlim(xlo, xhi)
    clipped = d3[(d3.mm_per_decade_lo < xlo) | (d3.mm_per_decade_hi > xhi)]
    for _, r_ in clipped.iterrows():
        side = xlo + 0.15 if r_.mm_per_decade_lo < xlo else xhi - 0.15
        ax.text(side, y[d3.index[d3.spp == r_.spp][0]], f"[{neg(f'{r_.mm_per_decade_lo:.1f}')}, {neg(f'{r_.mm_per_decade_hi:.1f}')}]",
                fontsize=5, va="center", ha="left" if r_.mm_per_decade_lo < xlo else "right", color="#555555")
    n_neg, n_excl_neg, n_excl_pos = int((d3.mm_per_decade < 0).sum()), int((d3.mm_per_decade_hi < 0).sum()), int((d3.mm_per_decade_lo > 0).sum())
    for k, v in [("n species", len(d3)), ("n negative", n_neg), ("n excl zero negative", n_excl_neg),
                 ("n excl zero positive", n_excl_pos), ("median mm/decade", float(d3.mm_per_decade.median()))]:
        note("fig-species-slopes", k, v)
    note("fig-species-slopes", "tier", src_label)
    ttl = (f"Species-specific wing-length trends under full controls (M3)\n{src_label}\n"
           f"{n_neg} of {len(d3)} slopes negative; {n_excl_neg} negative and {n_excl_pos} positive intervals exclude zero; "
           f"median {neg(f'{d3.mm_per_decade.median():.2f}')} mm/decade")
    if m3:
        ttl += f"\nblue line/band: pooled M3 year effect {neg(f'{m3['mmd']:.2f}')} [{neg(f'{m3['mmd_lo']:.2f}')}, {neg(f'{m3['mmd_hi']:.2f}')}] mm/decade"
    ax.set_title(ttl, fontsize=8.2, loc="left")
    handles = [Line2D([0], [0], marker="o", color=DEC, lw=0, label="negative slope"),
               Line2D([0], [0], marker="o", color=INC, lw=0, label="positive slope"),
               Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="black", label="interval includes zero"),
               Line2D([0], [0], marker="o", color="none", markerfacecolor="none", markeredgecolor="#8A8A8A", markersize=4, label="M0 baseline slope")]
    if d3_cmp is not None and len(d3_cmp):
        handles.append(Line2D([0], [0], marker="x", color="#2E5090", lw=0, alpha=0.55, label="lme4 Tier-1 M3 (comparison) slope"))
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig-species-slopes.png", bbox_inches="tight"); plt.close(fig)
    print(f"fig-species-slopes: {len(d3)} species, {n_neg} negative, {n_excl_neg}/{n_excl_pos} excl. zero (neg/pos)")
else:
    print("SKIPPED fig-species-slopes (fig_species_slopes.csv missing)")

# ---------------------------------------------------------------------------
# FIG 3 - per-species wing & bill trajectories, early vs late  [unchanged]
# ---------------------------------------------------------------------------
def per_species(metric):
    g = q.dropna(subset=[metric]).groupby(["Binomial", "period"])[metric].mean()
    return g.unstack().dropna()
fig, axes = plt.subplots(1, 2, figsize=(8, 4.3))
for ax, (metric, lab) in zip(axes, [("cwl", "Wing length"), ("Bill_width.mm.", "Bill width")]):
    wmat = per_species(metric)
    for sp, row in wmat.iterrows():
        e, l = row["1995–2006"], row["2013–2018"]
        ax.plot([0, 1], [e, l], color=(DEC if l < e else INC), alpha=0.6, lw=1)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["1995–2006", "2013–2018"])
    ax.set_ylabel(f"{lab} (mm)"); ax.set_title(lab, fontsize=10)
fig.legend([Line2D([0], [0], color=DEC), Line2D([0], [0], color=INC)],
           ["Decrease", "Increase"], loc="lower center", ncol=2, frameon=False,
           bbox_to_anchor=(0.5, -0.04))
fig.tight_layout(); fig.savefig(f"{OUT}/fig-trends.png", bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------------------
# FIG 4 - lnCVR forest plots, TWO panels (wing & bill); corrected SE  [unchanged]
# ---------------------------------------------------------------------------
def lncvr_table(metric):
    g = q.dropna(subset=[metric]).groupby(["Binomial", "period"])[metric].agg(m="mean", s="std", n="count").reset_index()
    piv = g.pivot(index="Binomial", columns="period").dropna()
    me, mc = piv[("m", "2013–2018")], piv[("m", "1995–2006")]
    se, sc = piv[("s", "2013–2018")], piv[("s", "1995–2006")]
    ne, nc = piv[("n", "2013–2018")], piv[("n", "1995–2006")]
    def corr(period):
        mm = np.log(g[g.period == period]["m"]); ss = np.log(g[g.period == period]["s"])
        ok = np.isfinite(mm) & np.isfinite(ss); return np.corrcoef(mm[ok], ss[ok])[0, 1]
    ce, cc = corr("2013–2018"), corr("1995–2006")
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

# ---------------------------------------------------------------------------
# SUPP FIG - wing records per year by contributor (stacked) + contributors per year
# ---------------------------------------------------------------------------
# Data: audit_sources.rds (P0): records_per_contributor_year_wing, contributor_table,
# contributors_per_year. The 12 largest contributors by wing records are named; the
# rest are pooled. Contributors whose wing records span both quartile periods
# (spans_both_wing) are marked with a dagger in the legend.
rcy = fd_csv("fig_records_per_contributor_year.csv")
ctab = fd_csv("fig_contributor_table.csv")
cpy = fd_csv("fig_contributors_per_year.csv")
top_k = 12
rank = ctab[ctab.n_wing > 0].sort_values("n_wing", ascending=False)
top = list(rank.Main_researcher.head(top_k))
n_other = int((rank.n_wing > 0).sum()) - top_k
rcy["contributor"] = np.where(rcy.Main_researcher.isin(top), rcy.Main_researcher, f"other ({n_other} contributors)")
piv = rcy.pivot_table(index="Year", columns="contributor", values="n_wing", aggfunc="sum", fill_value=0)
piv = piv.reindex(columns=top + [f"other ({n_other} contributors)"], fill_value=0)
years = np.arange(1995, 2019); piv = piv.reindex(years, fill_value=0)
spans = set(ctab.loc[ctab.spans_both_wing.astype(str).str.upper() == "TRUE", "Main_researcher"])
palette = list(plt.get_cmap("tab20").colors)[:top_k] + ["#BDBDBD"]
fig = plt.figure(figsize=(9, 6.2))
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.12)
ax = fig.add_subplot(gs[0]); bottom = np.zeros(len(years))
for c, colr in zip(piv.columns, palette):
    ax.bar(years, piv[c].to_numpy(), bottom=bottom, width=0.85, color=colr, edgecolor="white", linewidth=0.3,
           label=(c + " †") if c in spans else c)
    bottom += piv[c].to_numpy()
ax.set_ylabel("Wing-length records"); ax.set_xticks(years[::2]); ax.tick_params(labelbottom=False)
n_span, n_contrib_w, rec_span = S("n_span_both_wing", int), S("n_wing_contributors", int), S("wing_records_from_spanners", int)
ax.set_title(f"Wing records per year by contributor (Main_researcher; N = {int(piv.to_numpy().sum()):,} records, {n_contrib_w} contributors)\n"
             f"† contributor with wing records in both 1995–{EARLY_MAX} and {LATE_MIN}–2018: {n_span} of {n_contrib_w} ({rec_span:,} records)",
             fontsize=9)
ax.legend(loc="upper left", fontsize=6.6, frameon=False, ncol=2, title="contributor", title_fontsize=7)
ax2 = fig.add_subplot(gs[1], sharex=ax)
ax2.bar(cpy.Year, cpy.n_contributors_wing, width=0.85, color="#6E6E6E")
ax2.set_ylabel("Contributors"); ax2.set_xlabel("Year"); ax2.set_xticks(years[::2])
fig.savefig(f"{OUT}/fig-s-records-by-source.png", bbox_inches="tight"); plt.close(fig)
note("fig-s-records-by-source", "wing records", int(piv.to_numpy().sum()))
note("fig-s-records-by-source", "contributors with wing records", n_contrib_w)
note("fig-s-records-by-source", "contributors spanning both periods", n_span)
print(f"fig-s-records-by-source: {int(piv.to_numpy().sum())} records, top {top_k} named, {n_other} pooled; spanning marked: {sorted(spans)}")

# ---------------------------------------------------------------------------
# SUPP FIG - wing column populated, by year (protocol proxy)
# ---------------------------------------------------------------------------
wc = fd_csv("fig_wingcol_by_year.csv")
wpiv = wc.pivot_table(index="Year", columns="wing_col", values="prop", aggfunc="sum", fill_value=0)
wpiv = wpiv.reindex(columns=["right", "left", "generic"], fill_value=0).reindex(years, fill_value=0)
nyr = wc.groupby("Year").n_year.first().reindex(years)
wcol = {"right": "#2C7BB6", "left": "#ABD9E9", "generic": "#FDAE61"}
wlab = {"right": "Wing_length_right", "left": "Wing_length_left", "generic": "Wing_length (side unspecified)"}
fig, ax = plt.subplots(figsize=(8.6, 4.6)); bottom = np.zeros(len(years))
for c in wpiv.columns:
    ax.bar(years, wpiv[c].to_numpy(), bottom=bottom, width=0.85, color=wcol[c], edgecolor="white", linewidth=0.3, label=wlab[c])
    bottom += wpiv[c].to_numpy()
for yr, n in nyr.items():
    if pd.notna(n):
        ax.text(yr, 1.01, f"{int(n)}", ha="center", va="bottom", fontsize=6, rotation=90)
ax.set_ylim(0, 1.16); ax.set_yticks([0, .25, .5, .75, 1]); ax.set_yticklabels(["0", "25", "50", "75", "100"])
ax.set_ylabel("Share of wing records (%)"); ax.set_xlabel("Year"); ax.set_xticks(years[::2])
p_right_early, p_generic_late = S("prop_right_early"), S("prop_generic_late")
ax.set_title(f"Which wing column supplied the coalesced wing length, by year (numbers above bars = records)\n"
             f"protocol proxy: {100 * p_right_early:.0f}% right-wing column in 1995–{EARLY_MAX} vs "
             f"{100 * p_generic_late:.0f}% side-unspecified column in {LATE_MIN}–2018", fontsize=9)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3, frameon=False, fontsize=8, title="column populated", title_fontsize=8)
fig.tight_layout(); fig.savefig(f"{OUT}/fig-s-wingcol.png", bbox_inches="tight"); plt.close(fig)
note("fig-s-wingcol", "prop right early", p_right_early); note("fig-s-wingcol", "prop generic late", p_generic_late)

pd.DataFrame(panel_numbers).to_csv(os.path.join(FD, "fig_panel_numbers.csv"), index=False)
print("figures written to", OUT, "| panel numbers ->", os.path.join(FD, "fig_panel_numbers.csv"))
