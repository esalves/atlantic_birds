import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec

OUT="Manuscript/images"
df = pd.read_csv("Analysis/ATLANTIC_BIRD_TRAITS_completed_2018_11_d05.csv", low_memory=False)
f = df[(df.Year>=1990)&(df.Age=="Adult")&(df.Sex!="Unknown")&(df.Order=="Passeriformes")&(df.AtlanticForests_20km_Buffer=="inside the 20 km polygon")].copy()
f["cwl"]=f["Wing_length_right.mm."].combine_first(f["Wing_length_left.mm."]).combine_first(f["Wing_length.mm."])
vc=f.groupby("Binomial").size(); f=f[f.Binomial.isin(vc[vc>=50].index)]
rng=f.groupby("Binomial").Year.agg(lambda s:s.max()-s.min()); f=f[f.Binomial.isin(rng[rng>=10].index)]
f["Binomial"]=f["Binomial"].str.replace(" ","_")
print("species",f.Binomial.nunique(),"records",len(f))

# periods (1st and 3rd quartile windows used in the Rmd)
q=f[(f.Year<=2006)|(f.Year>=2013)].copy()
q["period"]=np.where(q.Year<=2006,"1990–2006","2013–2018")

plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":150})
DEC="#CC6666"; INC="#5B7FB5"

# ---------- FIG 1: map by period ----------
fig=plt.figure(figsize=(8,4.2))
gs=gridspec.GridSpec(1,2,wspace=0.18)
for i,p in enumerate(["1990–2006","2013–2018"]):
    ax=fig.add_subplot(gs[0,i]); d=q[(q.period==p)].dropna(subset=["cwl"])
    sc=ax.scatter(d.Longitude_decimal_degrees,d.Latitude_decimal_degrees,c=d.cwl,
                  cmap="plasma",s=9,alpha=0.55,edgecolors="none",vmin=q.cwl.quantile(.02),vmax=q.cwl.quantile(.98))
    ax.set_title(p,fontsize=10); ax.set_xlabel("Longitude (°)")
    if i==0: ax.set_ylabel("Latitude (°)")
    ax.set_aspect("equal","box")
cb=fig.colorbar(sc,ax=fig.axes,fraction=0.025,pad=0.02); cb.set_label("Wing length (mm)")
fig.suptitle("Sampling localities of 68 Atlantic Forest passerine species",fontsize=11,y=0.99)
fig.savefig(f"{OUT}/fig-map.png",bbox_inches="tight"); plt.close(fig)

# ---------- species means per period ----------
def per_species(metric):
    g=q.dropna(subset=[metric]).groupby(["Binomial","period"])[metric].agg(["mean","std","count"]).reset_index()
    w=g.pivot(index="Binomial",columns="period",values="mean")
    w=w.dropna()
    return g,w
gwing,wwing=per_species("cwl")
gbill,wbill=per_species("Bill_width.mm.")

# ---------- FIG 2: trends ----------
fig,axes=plt.subplots(1,2,figsize=(8,4.3))
for ax,(w,lab,unit) in zip(axes,[(wwing,"Wing length","mm"),(wbill,"Bill width","mm")]):
    dec=inc=0
    for sp,row in w.iterrows():
        e,l=row["1990–2006"],row["2013–2018"]; col=DEC if l<e else INC
        dec+=l<e; inc+=l>=e
        ax.plot([0,1],[e,l],color=col,alpha=0.6,lw=1)
    ax.set_xticks([0,1]); ax.set_xticklabels(["1990–2006","2013–2018"])
    ax.set_ylabel(f"{lab} ({unit})"); ax.set_title(lab,fontsize=10)
    print(lab,"dec",dec,"inc",inc)
from matplotlib.lines import Line2D
fig.legend([Line2D([0],[0],color=DEC),Line2D([0],[0],color=INC)],["Decrease","Increase"],
           loc="lower center",ncol=2,frameon=False,bbox_to_anchor=(0.5,-0.04))
fig.tight_layout(); fig.savefig(f"{OUT}/fig-trends.png",bbox_inches="tight"); plt.close(fig)

# ---------- lnCVR (corrected CI uses sqrt of sampling variance) ----------
def lncvr_table(metric):
    g=q.dropna(subset=[metric]).groupby(["Binomial","period"])[metric].agg(m="mean",s="std",n="count").reset_index()
    piv=g.pivot(index="Binomial",columns="period")
    piv=piv.dropna()
    me=piv[("m","2013–2018")]; mc=piv[("m","1990–2006")]
    se=piv[("s","2013–2018")]; sc=piv[("s","1990–2006")]
    ne=piv[("n","2013–2018")]; nc=piv[("n","1990–2006")]
    # correlation across species of log mean vs log sd, per period
    def corr(period):
        mm=np.log(g[g.period==period]["m"]); ss=np.log(g[g.period==period]["s"])
        ok=np.isfinite(mm)&np.isfinite(ss); return np.corrcoef(mm[ok],ss[ok])[0,1]
    cor_e=corr("2013–2018"); cor_c=corr("1990–2006")
    lncvr=np.log((se/me)/(sc/mc)) + 1/(2*(ne-1)) - 1/(2*(nc-1))
    def s2(mc,sc,nc,corc,me,se,ne,core):
        ac=sc**2/(nc*mc**2); bc=1/(2*(nc-1)); cc=2*corc*np.sqrt(ac*bc)
        ae=se**2/(ne*me**2); be=1/(2*(ne-1)); ce=2*core*np.sqrt(ae*be)
        return ac+bc-cc+ae+be-ce
    var=s2(mc,sc,nc,cor_c,me,se,ne,cor_e)
    t=pd.DataFrame({"sp":piv.index,"lncvr":lncvr.values,"var":var.values}).dropna()
    t["se"]=np.sqrt(t["var"])  # CORRECTED: standard error, not the variance
    t=t.sort_values("lncvr")
    return t
tw=lncvr_table("cwl")
# random-effects (DL) mean for annotation
def re_mean(t):
    w=1/t["var"]; mu=(w*t["lncvr"]).sum()/w.sum()
    Q=(w*(t["lncvr"]-mu)**2).sum(); k=len(t); C=w.sum()-(w**2).sum()/w.sum()
    tau2=max(0,(Q-(k-1))/C); wr=1/(t["var"]+tau2)
    mur=(wr*t["lncvr"]).sum()/wr.sum(); se=np.sqrt(1/wr.sum())
    return mur,mur-1.96*se,mur+1.96*se,k
mu,lo,hi,k=re_mean(tw); print(f"wing lnCVR RE mean {mu:.3f} [{lo:.3f},{hi:.3f}] k={k}")

# ---------- FIG 3: forest plot wing lnCVR ----------
fig,ax=plt.subplots(figsize=(5.2,8))
y=np.arange(len(tw))
ax.errorbar(tw["lncvr"],y,xerr=1.96*tw["se"],fmt="o",ms=3,lw=0.7,color="#333333",ecolor="#999999",capsize=1.5)
ax.axvline(0,ls="--",color="grey",lw=0.8)
ax.axvspan(lo,hi,color="#5B7FB5",alpha=0.18)
ax.axvline(mu,color="#2E5090",lw=1.4)
ax.set_yticks(y); ax.set_yticklabels([s.replace("_"," ") for s in tw["sp"]],fontsize=5.5,style="italic")
ax.set_xlabel("lnCVR (wing length), late vs early  ·  95% CI")
ax.set_title(f"Among-individual variability in wing length\nmeta-analytic mean = {mu:.2f} [{lo:.2f}, {hi:.2f}]",fontsize=9)
fig.tight_layout(); fig.savefig(f"{OUT}/fig-variability.png",bbox_inches="tight"); plt.close(fig)
print("figures written")
