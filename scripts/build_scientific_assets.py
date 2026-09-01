"""Build the Cinch methods-paper figures and frozen SPN534 presentation assets.

This script is deliberately presentation-only.  It reads frozen MI/NMI/Neff and
legacy distance outputs; it never recomputes or overwrites association scores.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import html
import json
import math
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


COLORS = {
    "navy": "#16324F",
    "blue": "#2C73B9",
    "cyan": "#4AA8B7",
    "green": "#3A9D72",
    "orange": "#E07A5F",
    "gold": "#E9B44C",
    "red": "#C94C4C",
    "purple": "#7353BA",
    "grey": "#79838D",
    "light": "#F4F7F9",
    "ink": "#1D2731",
}
CHANNEL_COLORS = {"PP": "#2676B8", "PT": "#E69F00", "TP": "#8B5CF6", "TT": "#D84A4A"}


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 130,
            "savefig.dpi": 220,
            "savefig.bbox": "tight",
            "axes.grid": True,
            "grid.alpha": 0.16,
        }
    )


def save(fig: plt.Figure, base: Path) -> None:
    base.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(base.with_suffix(f".{ext}"), facecolor="white")
    plt.close(fig)


def box(ax, xy, width, height, text, fc="white", ec=None, fontsize=10, lw=1.5):
    ec = ec or COLORS["navy"]
    p = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.015,rounding_size=0.02",
        facecolor=fc,
        edgecolor=ec,
        linewidth=lw,
    )
    ax.add_patch(p)
    ax.text(xy[0] + width / 2, xy[1] + height / 2, text, ha="center", va="center", fontsize=fontsize)
    return p


def arrow(ax, start, end, color=None, rad=0.0):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=1.6,
            color=color or COLORS["navy"],
            connectionstyle=f"arc3,rad={rad}",
        )
    )


def conceptual_figures(out: Path) -> None:
    # FIG01 — why Cinch exists
    fig, ax = plt.subplots(figsize=(11, 4.8)); ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    box(ax, (0.03, .57), .25, .25, "Gene presence/absence\nsees locus carriage", "#E7F1FA", COLORS["blue"], 12)
    box(ax, (.375, .57), .25, .25, "Core loci can be\npresent in every genome", "#F5EFE3", COLORS["gold"], 12)
    box(ax, (.72, .57), .25, .25, "Allele/type states reveal\nhidden dependence", "#FBE8E7", COLORS["red"], 12)
    arrow(ax, (.28, .695), (.375, .695)); arrow(ax, (.625, .695), (.72, .695))
    ax.text(.5, .35, "Cinch tests information dependence across presence and type states", ha="center", fontsize=15, weight="bold", color=COLORS["navy"])
    ax.text(.5, .22, "Then asks whether the signal is long-range and recurrent across HC69 backgrounds.", ha="center", fontsize=11)
    ax.set_title("FIG01  Why Cinch exists", loc="left", weight="bold", fontsize=17)
    save(fig, out / "FIG01_WHY_CINCH")

    # FIG02 — state model
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    samples = np.arange(12)
    presence = np.array([1,1,0,1,1,1,0,1,1,1,1,0])
    types = np.array([1,2,0,1,3,2,0,3,1,2,3,0])
    axes[0].imshow(presence[None, :], aspect="auto", cmap=LinearSegmentedColormap.from_list("p", ["#EEEEEE", COLORS["blue"]]), vmin=0, vmax=1)
    axes[0].set_title("Presence state  P ∈ {0,1}"); axes[0].set_yticks([0], ["locus A"]); axes[0].set_xticks(samples, [f"S{i+1}" for i in samples], rotation=60)
    axes[1].imshow(types[None, :], aspect="auto", cmap="Set2", vmin=0, vmax=3)
    axes[1].set_title("Type state  T ∈ {a₁,a₂,…}"); axes[1].set_yticks([0], ["locus A"]); axes[1].set_xticks(samples, [f"S{i+1}" for i in samples], rotation=60)
    fig.suptitle("FIG02  One biological-state representation, two resolutions", fontsize=17, weight="bold")
    fig.text(.5, .03, "Missing/absent and nominal allele categories remain explicit; allele IDs are categories, not numerical distances.", ha="center")
    fig.tight_layout(rect=(0,.07,1,.9)); save(fig, out / "FIG02_STATE_MODEL")

    # FIG03 — three association scales / four channels
    fig, ax = plt.subplots(figsize=(11, 5)); ax.set_axis_off(); ax.set_xlim(0,1); ax.set_ylim(0,1)
    coords = {"PP": (.12,.64), "PT": (.38,.64), "TP": (.64,.64), "TT": (.86,.64)}
    for ch,(x,y) in coords.items():
        fc = CHANNEL_COLORS[ch] + "22"
        box(ax, (x-.09,y-.10), .18,.20, ch, fc, CHANNEL_COLORS[ch], 16, 2)
    ax.text(.5,.88,"Four state channels",ha="center",fontsize=16,weight="bold")
    box(ax,(.12,.20),.22,.18,"MIraw / NMI\nassociation strength","#E7F1FA",COLORS["blue"],12)
    box(ax,(.39,.20),.22,.18,"HC69 recurrence\nbackground support","#EDF7F1",COLORS["green"],12)
    box(ax,(.66,.20),.22,.18,"Order distance\nlocal-linkage axis","#FBEFE6",COLORS["orange"],12)
    for x in (.23,.50,.77): arrow(ax,(x,.57),(x,.39))
    ax.set_title("FIG03  Cinch separates association, population recurrence and linkage",loc="left",fontsize=17,weight="bold")
    save(fig,out/"FIG03_THREE_SCALES")

    # FIG04 — MI and driver cells
    obs=np.array([[38,2,1],[3,29,4],[1,5,27]],float); exp=np.outer(obs.sum(1),obs.sum(0))/obs.sum(); res=(obs-exp)/np.sqrt(exp)
    fig,axes=plt.subplots(1,3,figsize=(12,4.2))
    im=axes[0].imshow(obs,cmap="Blues"); axes[0].set_title("Observed state pairs")
    axes[1].imshow(exp,cmap="Greys"); axes[1].set_title("Expected under independence")
    im2=axes[2].imshow(res,cmap="RdBu_r",norm=TwoSlopeNorm(vcenter=0)); axes[2].set_title("Pearson residual = driver map")
    for ax in axes:
        ax.set_xticks(range(3),["B1","B2","B3"]); ax.set_yticks(range(3),["A1","A2","A3"])
    fig.colorbar(im2,ax=axes[2],fraction=.046); fig.suptitle("FIG04  MI measures total dependence; cell residuals identify the states that drive it",fontsize=15,weight="bold")
    fig.tight_layout(rect=(0,0,1,.9)); save(fig,out/"FIG04_MI_AND_DRIVERS")

    # FIG05 — full workflow
    fig,ax=plt.subplots(figsize=(13,4.7)); ax.set_axis_off(); ax.set_xlim(0,1); ax.set_ylim(0,1)
    labels=[("Assemblies + GFF\n+ core profiles",COLORS["grey"]),("P/T state\nmatrices",COLORS["blue"]),("PP · PT · TP · TT\nMIraw / NMI",COLORS["purple"]),("Order-distance\nlong-range filter",COLORS["orange"]),("HC69 cross-background\nrecurrence + Neff",COLORS["green"]),("ARACNE\nredundancy pruning",COLORS["red"]),("Frozen network\n223 edges",COLORS["navy"])]
    xs=np.linspace(.02,.86,len(labels))
    for i,((lab,c),x) in enumerate(zip(labels,xs)):
        box(ax,(x,.40),.12,.24,lab,c+"20",c,9.5,1.8)
        if i<len(labels)-1: arrow(ax,(x+.12,.52),(xs[i+1],.52),c)
    ax.text(.5,.20,"HC69 is a recurrence filter, not a corrected p/q-value.  ARACNE acts only after filtering.",ha="center",fontsize=11,weight="bold")
    ax.set_title("FIG05  Frozen Cinch WGS workflow",loc="left",fontsize=17,weight="bold")
    save(fig,out/"FIG05_WGS_PIPELINE")


def distance_figure(source: Path, repo: Path, out: Path) -> pd.DataFrame:
    p = source / "track_B_retrospective/figures/data/original_distance_diagnostic/SPN534_ORIGINAL_DISTANCE_DIAGNOSTIC_ALL_PAIRS.parquet"
    df = pd.read_parquet(p)
    stats=[]
    total=len(df)
    both=df["order_dist"].notna() & df["bp_dist"].notna()
    for key,col,unit in [("phylo","phylo_dist","legacy tree units"),("order","order_dist","genes"),("bp","bp_dist","bp")]:
        m=df[col].notna(); rho,pv=spearmanr(df.loc[m,col],df.loc[m,"MIraw"])
        x=df.loc[m,col].to_numpy(float); y=df.loc[m,"MIraw"].to_numpy(float)
        # Natural-width summaries; plotting uses log-spaced bins for genomic distance.
        if key=="phylo": edges=np.linspace(np.nanmin(x),np.nanmax(x),81)
        else: edges=np.unique(np.r_[0,np.geomspace(max(1,np.nanmin(x[x>0]) if np.any(x>0) else 1),np.nanmax(x)+1,81)])
        b=np.digitize(x,edges)-1; rows=[]
        for bi in range(len(edges)-1):
            yy=y[b==bi]; xx=x[b==bi]
            if len(yy): rows.append({"distance_type":key,"bin":bi,"N":len(yy),"distance_median":np.median(xx),"MI_median":np.median(yy),"MI_mean":np.mean(yy),"MI_Q25":np.quantile(yy,.25),"MI_Q75":np.quantile(yy,.75),"MI_Q95":np.quantile(yy,.95),"MI_Q99":np.quantile(yy,.99),"MI_Q99_9":np.quantile(yy,.999)})
        pd.DataFrame(rows).to_csv(repo/"examples/spn534/plot_data"/f"DISTANCE_BINS_{key.upper()}.tsv",sep="\t",index=False)
        stats.append({"distance":key,"definition":unit,"available_pairs":int(m.sum()),"missing_pairs":int((~m).sum()),"availability_fraction":float(m.mean()),"spearman_rho_MIraw":rho,"spearman_p_descriptive":pv})
    od=df.loc[both,"order_dist"].to_numpy(float); bp=df.loc[both,"bp_dist"].to_numpy(float)
    rho_ob,p_ob=spearmanr(od,bp)
    table=pd.DataFrame(stats)
    table["total_pairs"]=total
    table["both_order_bp_available"]=int(both.sum())
    table["order_only_available"]=int((df.order_dist.notna() & df.bp_dist.isna()).sum())
    table["bp_only_available"]=int((df.bp_dist.notna() & df.order_dist.isna()).sum())
    table["neither_order_nor_bp_available"]=int((df.bp_dist.isna() & df.order_dist.isna()).sum())
    table["order_bp_spearman_rho_when_both"]=rho_ob
    table["order_bp_spearman_p_descriptive"]=p_ob
    table["different_contig_or_no_reliable_coobservation_loss"]=int((df.n_same_contig_observations.fillna(0)==0).sum())
    table["insufficient_same_contig_observations_1_to_4"]=int(df.n_same_contig_observations.fillna(0).between(1,4).sum())
    diag=repo/"examples/spn534/diagnostics"; diag.mkdir(parents=True,exist_ok=True)
    table.to_csv(diag/"TABLE_DISTANCE_COMPARISON.tsv",sep="\t",index=False)

    fig,axes=plt.subplots(2,3,figsize=(14,8),gridspec_kw={"height_ratios":[2.2,1]})
    settings=[("phylo_dist","Phylogenetic distance","legacy tree units",False),("order_dist","Gene-order distance","genes",True),("bp_dist","Physical distance","bp",True)]
    for j,(col,title,unit,use_log) in enumerate(settings):
        m=df[col].notna(); x=df.loc[m,col].to_numpy(float); y=df.loc[m,"MIraw"].to_numpy(float); xx=np.log10(x+1) if use_log else x
        hb=axes[0,j].hexbin(xx,y,gridsize=75,mincnt=1,bins="log",cmap="mako" if "mako" in plt.colormaps() else "viridis",linewidths=0)
        rho=table.loc[table.distance==col.split('_')[0],"spearman_rho_MIraw"].iloc[0]
        axes[0,j].set(title=title,xlabel=(f"log10({unit} + 1)" if use_log else unit),ylabel="MIraw (nats)")
        axes[0,j].text(.98,.96,f"N = {m.sum():,}\nSpearman ρ = {rho:.3f}",transform=axes[0,j].transAxes,ha="right",va="top",bbox=dict(fc="white",alpha=.85,ec="#BBBBBB"))
        fig.colorbar(hb,ax=axes[0,j],label="log10 pair density")
    vals=[df.phylo_dist.notna().sum(),df.order_dist.notna().sum(),df.bp_dist.notna().sum()]
    axes[1,0].bar(["phylo","order","bp"],np.array(vals)/total,color=[COLORS["blue"],COLORS["green"],COLORS["orange"]]); axes[1,0].set_ylim(0,1); axes[1,0].set_ylabel("fraction of all pairs"); axes[1,0].set_title("Coordinate availability")
    axes[1,1].hexbin(np.log10(od+1),np.log10(bp+1),gridsize=65,mincnt=1,bins="log",cmap="viridis",linewidths=0); axes[1,1].set(xlabel="log10(order + 1)",ylabel="log10(bp + 1)",title=f"Order–bp agreement (ρ={rho_ob:.3f})")
    axes[1,2].axis("off"); axes[1,2].text(0,.9,"Why order distance is retained",fontsize=13,weight="bold",color=COLORS["navy"]); axes[1,2].text(0,.70,"• same exact coordinate evidence as bp distance",fontsize=11); axes[1,2].text(0,.55,"• invariant to indel-expanded intergenic length",fontsize=11); axes[1,2].text(0,.40,"• transferable across assemblies",fontsize=11); axes[1,2].text(0,.20,"Different-contig pairs remain NA — never ‘very far’.",fontsize=10,color=COLORS["red"])
    fig.suptitle("FIG06  Legacy phylogenetic, gene-order and bp distance diagnostics",fontsize=17,weight="bold")
    fig.tight_layout(rect=(0,0,1,.96)); save(fig,out/"FIG06_DISTANCE_COMPARISON")
    return table


def data_figures(repo: Path, source: Path, out: Path) -> None:
    plot=repo/"examples/spn534/plot_data"; results=repo/"examples/spn534/results"
    # The original fine-resolution fixed-width-bin mean-MI plot is retained as
    # the lead distance visual; FIG06_DISTANCE_COMPARISON adds density and
    # availability diagnostics rather than replacing it.
    for ext in ("png","pdf","svg"):
        shutil.copy2(
            source/f"track_B_retrospective/figures/FIG_SPN534_ORIGINAL_MIRAW_DISTANCE_MEAN_SCATTER_1X3.{ext}",
            out/f"FIG06_DISTANCE_MEAN_CURVES.{ext}",
        )
    # FIG07 — order background envelopes plus threshold
    fig,axes=plt.subplots(2,2,figsize=(11,8),sharex=True)
    for ax,ch in zip(axes.flat,["PP","PT","TP","TT"]):
        d=pd.read_csv(plot/f"ORDER_BACKGROUND_{ch}.tsv",sep="\t")
        x=np.log10(d.distance_median+1)
        ax.plot(x,d.median_MIraw,color=CHANNEL_COLORS[ch],lw=2,label="median")
        ax.fill_between(x,d.Q25,d.Q75,color=CHANNEL_COLORS[ch],alpha=.18,label="IQR")
        ax.plot(x,d.Q99_9,color=COLORS["ink"],lw=1.2,ls="--",label="Q99.9")
        ax.axvline(np.log10(101),color=COLORS["red"],lw=1.4,ls=":")
        ax.set(title=ch,ylabel="MIraw (nats)"); ax.legend(fontsize=8)
    for ax in axes[-1]: ax.set_xlabel("log10(order distance + 1), genes")
    fig.suptitle("FIG07  Channel-specific order-distance backgrounds and frozen 100-gene threshold",fontsize=16,weight="bold")
    fig.tight_layout(rect=(0,0,1,.95)); save(fig,out/"FIG07_ORDER_BACKGROUND")

    # FIG08 — HC hierarchy selection
    h=pd.read_csv(plot/"HC_LEVEL_SCAN.tsv",sep="\t",na_values="NA")
    fig,axes=plt.subplots(2,2,figsize=(11,7))
    axes[0,0].plot(h.HC_level,h.n_clusters,color=COLORS["blue"]); axes[0,0].set(ylabel="clusters",title="Resolution")
    axes[0,1].plot(h.HC_level,h.silhouette,color=COLORS["green"]); axes[0,1].set(ylabel="silhouette",title="Separation")
    axes[1,0].plot(h.HC_level,h.singleton_fraction,color=COLORS["orange"]); axes[1,0].set(ylabel="singleton fraction",title="Fragmentation")
    axes[1,1].plot(h.HC_level,h.NMI_to_previous,color=COLORS["purple"],label="NMI"); axes[1,1].plot(h.HC_level,h.ARI_to_previous,color=COLORS["red"],label="ARI"); axes[1,1].legend(); axes[1,1].set(ylabel="agreement",title="Neighbour-level stability")
    for ax in axes.flat: ax.axvline(69,color=COLORS["red"],ls="--",lw=1.5); ax.set_xlabel("HC threshold (allele differences)")
    fig.suptitle("FIG08  HC69 selected as a population-background resolution",fontsize=16,weight="bold")
    fig.tight_layout(rect=(0,0,1,.95)); save(fig,out/"FIG08_HC_SELECTION")

    # FIG09 — Neff derivation with exact examples
    fig,axes=plt.subplots(1,3,figsize=(13,4.2))
    examples=[("single background",[10]),("gyrB–aguA",[1,3,3]),("distributed support",[10,9,8,7,6])]
    for ax,(name,counts) in zip(axes,examples):
        c=np.asarray(counts,float); p=c/c.sum(); neff=1/np.sum(p*p)
        ax.bar(range(1,len(c)+1),c,color=COLORS["green"]); ax.set(title=f"{name}\nNeff = {neff:.2f}",xlabel="HC69 background",ylabel="driver support")
    fig.suptitle("FIG09  Effective number of supporting backgrounds",fontsize=16,weight="bold")
    fig.text(.5,.015,"pₕ = nₕ / Σₕnₕ    and    Neff = 1 / Σₕpₕ².  Neff quantifies dispersion; it is not a sample-size correction.",ha="center",fontsize=11)
    fig.tight_layout(rect=(0,.06,1,.92)); save(fig,out/"FIG09_NEFF_DERIVATION")

    # FIG10 — development contrast: rsuA_2–glyS (single HC) vs frozen gyrB–aguA
    rsu=pd.read_csv(source/"raw_type_channels/figures/data/tt_contingency_heatmaps/01_rsuA_2__glyS.tsv",sep="\t")
    mat1=rsu.pivot(index="state_A",columns="state_B",values="pearson_residual")
    npz=np.load(source/"raw_type_channels/matrices/SPN534_UNIFIED_STATE_MATRIX.npz",allow_pickle=True)
    loci=[str(x) for x in npz["loci"]]; ts=npz["type_state"]
    ia,ib=loci.index("gyrB"),loci.index("aguA"); a=ts[:,ia]; b=ts[:,ib]; valid=(a>0)&(b>0)
    ta=pd.crosstab(pd.Series(a[valid],name="A"),pd.Series(b[valid],name="B")); obs=ta.to_numpy(float); exp=np.outer(obs.sum(1),obs.sum(0))/obs.sum(); resid=(obs-exp)/np.sqrt(np.where(exp>0,exp,np.nan)); mat2=pd.DataFrame(resid,index=[f"type_{i}" for i in ta.index],columns=[f"type_{i}" for i in ta.columns])
    cdir=plot/"contingency"; cdir.mkdir(exist_ok=True); mat2.to_csv(cdir/"GYRB_AGUA_PEARSON_RESIDUAL.tsv",sep="\t")
    fig,axes=plt.subplots(1,2,figsize=(13,5.5))
    vmax=max(np.nanpercentile(np.abs(mat1.to_numpy()),99),np.nanpercentile(np.abs(mat2.to_numpy()),99)); norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax)
    for ax,mat,title in [(axes[0],mat1,"rsuA_2–glyS\nTT-high but dominant driver in one HC69; Neff = 1"),(axes[1],mat2,"gyrB–aguA\nfrozen TT edge; HC support 1:3:3; Neff = 2.58")]:
        im=ax.imshow(mat,cmap="RdBu_r",norm=norm,aspect="auto"); ax.set_title(title); ax.set_xlabel(mat.columns.name or "locus B state"); ax.set_ylabel(mat.index.name or "locus A state");
        sx=max(1,len(mat.columns)//12); sy=max(1,len(mat.index)//12)
        ax.set_xticks(range(0,len(mat.columns),sx),mat.columns[::sx],rotation=90,fontsize=6)
        ax.set_yticks(range(0,len(mat.index),sy),mat.index[::sy],fontsize=6)
        if "type_46" in mat.index and "type_51" in mat.columns:
            yi=list(mat.index).index("type_46"); xi=list(mat.columns).index("type_51")
            ax.add_patch(plt.Rectangle((xi-.5,yi-.5),1,1,fill=False,edgecolor="#FFD23F",linewidth=2.5))
            ax.annotate("driver",(xi,yi),xytext=(xi-8,yi-7),arrowprops=dict(arrowstyle="->",color="#8A5A00"),fontsize=8,color="#8A5A00")
    fig.colorbar(im,ax=axes.ravel().tolist(),label="Pearson residual",fraction=.025)
    fig.suptitle("FIG10  High TT information is not sufficient: background recurrence changes interpretation",fontsize=15,weight="bold")
    fig.subplots_adjust(top=.82,bottom=.24,wspace=.25); save(fig,out/"FIG10_TT_BACKGROUND_CONTRAST")

    # FIG11 — exact all-pair Neff/order diagnostic, copied under paper numbering
    for ext in ("png","pdf","svg"):
        shutil.copy2(repo/f"examples/spn534/figures/NEFF_VS_ORDER.{ext}",out/f"FIG11_NEFF_VS_ORDER.{ext}")

    # FIG12 — Coinfinder concordance
    sm=pd.read_csv(source/"raw_type_channels/results/TABLE_SPNA_HIGH_PAIR_COINFINDER_SUMMARY.tsv",sep="\t")
    d=dict(zip(sm.metric,sm.value)); n=int(d["candidate_pairs"]); direct=int(d["Coinfinder_association_supported"]+d["Coinfinder_dissociation_supported"]); comp=int(d["same_published_component"])
    fig,axes=plt.subplots(1,2,figsize=(10,4.5))
    axes[0].bar(["Direct published\nCoinfinder edge","Not direct"],[direct,n-direct],color=[COLORS["blue"],"#D9DEE3"]); axes[0].set_ylabel("high-MI candidate pairs"); axes[0].set_title(f"Direct edge concordance\n{direct}/{n} ({direct/n:.1%})")
    axes[1].bar(["Same published\ncomponent","Outside component"],[comp,n-comp],color=[COLORS["green"],"#D9DEE3"]); axes[1].set_title(f"Component concordance\n{comp}/{n} ({comp/n:.1%})")
    fig.suptitle("FIG12  Published Coinfinder comparison is a concordance benchmark, not causal truth",fontsize=15,weight="bold")
    fig.tight_layout(rect=(0,0,1,.9)); save(fig,out/"FIG12_COINFINDER_CONCORDANCE")

    # FIG13 and FIG14 — exact existing frozen-data graphics
    for ext in ("png","pdf","svg"):
        shutil.copy2(source/f"raw_type_channels/figures/FIG_SPN_PP_VS_TT_RAW_MI.{ext}",out/f"FIG13_PP_VS_TT.{ext}")
        shutil.copy2(source/f"final_neff_filtered/figures/FIG_NEFF_FILTER03_PP_TT_TREE_PHANDANGO.{ext}",out/f"FIG14_TREE_HC_STATE.{ext}")

    # FIG15 — frozen attrition
    for ext in ("png","pdf","svg"):
        shutil.copy2(repo/f"examples/spn534/figures/FILTER_ATTRITION.{ext}",out/f"FIG15_FILTER_ATTRITION.{ext}")


def build_network(repo: Path, out: Path) -> tuple[nx.Graph, pd.DataFrame, pd.DataFrame]:
    final=pd.read_csv(repo/"examples/spn534/results/FINAL_CANDIDATES.tsv",sep="\t",na_values="NA")
    final["pair_class"]=final.channel.replace({"PT":"P↔T","TP":"P↔T"})
    G=nx.Graph()
    for r in final.itertuples(index=False):
        attrs={"edge_id":str(r.edge_id),"channel":str(r.channel),"pair_class":str(r.pair_class),"orientation":str(r.orientation),"MIraw":float(r.MIraw),"NMI":float(r.NMI),"order_distance":float(r.order_distance),"HC_supported":int(r.HC_supported),"Neff":float(r.Neff),"driver_A":str(r.enriched_driver_A_state),"driver_B":str(r.enriched_driver_B_state),"annotation_A":str(r.annotation_A),"annotation_B":str(r.annotation_B)}
        G.add_edge(str(r.locus_A),str(r.locus_B),**attrs)
        for locus,ann in [(str(r.locus_A),str(r.annotation_A)),(str(r.locus_B),str(r.annotation_B))]:
            if locus not in G.nodes: G.add_node(locus,annotation=ann)
            elif G.nodes[locus].get("annotation") in (None,"nan","NA"): G.nodes[locus]["annotation"]=ann
    communities=list(nx.community.greedy_modularity_communities(G,weight="NMI")); cmap={n:i+1 for i,c in enumerate(communities) for n in c}
    # Lay out the large component for readable structure and pack tiny components
    # along the lower margin, so satellites do not crush the main graph.
    components=sorted(nx.connected_components(G),key=len,reverse=True); pos={}
    main=G.subgraph(components[0]); p0=nx.spring_layout(main,seed=42,weight="NMI",iterations=350,k=.28)
    for n,(x,y) in p0.items(): pos[n]=(x,0.90*y+0.18)
    satellites=components[1:]
    for ci,c in enumerate(satellites):
        sg=G.subgraph(c); local=nx.circular_layout(sg); cx=-.84+1.68*(ci/max(1,len(satellites)-1))
        for n,(x,y) in local.items(): pos[n]=(cx+.055*x,-1.08+.055*y)
    for n in G:
        G.nodes[n].update(degree=int(G.degree(n)),community=int(cmap[n]),x=float(pos[n][0]),y=float(pos[n][1]))
    nd=pd.DataFrame([{"locus":n,**G.nodes[n]} for n in G.nodes]).sort_values(["community","degree"],ascending=[True,False])
    ed=pd.DataFrame([{"locus_A":u,"locus_B":v,**a} for u,v,a in G.edges(data=True)])
    net=repo/"examples/spn534/network"; net.mkdir(parents=True,exist_ok=True)
    nd.to_csv(net/"NETWORK_NODES.tsv",sep="\t",index=False); ed.to_csv(net/"NETWORK_EDGES.tsv",sep="\t",index=False)
    nx.write_graphml(G,net/"CINCH.graphml"); nx.write_gexf(G,net/"CINCH.gexf")
    # FIG16 network
    fig,ax=plt.subplots(figsize=(12,10)); ax.set_axis_off()
    node_colors=[plt.cm.tab20((cmap[n]-1)%20) for n in G]; sizes=[18+9*math.sqrt(G.degree(n)) for n in G]
    nx.draw_networkx_edges(G,pos,ax=ax,width=[.25+1.8*G.edges[e]["NMI"] for e in G.edges],alpha=.27,edge_color=[CHANNEL_COLORS[G.edges[e]["channel"]] for e in G.edges])
    nx.draw_networkx_nodes(G,pos,ax=ax,node_size=sizes,node_color=node_colors,edgecolors="white",linewidths=.25)
    hubs=sorted(G.degree,key=lambda x:x[1],reverse=True)[:6]
    for li,(n,_) in enumerate(hubs):
        label=G.nodes[n].get("annotation")
        if label in (None,"nan","Uncharacterised protein"): label=n
        x,y=pos[n]; dx=.13 if x<0 else -.13; dy=.10*((li%3)-1)
        ax.annotate(str(label),xy=(x,y),xytext=(x+dx,y+dy),ha="left" if dx>0 else "right",fontsize=7,arrowprops=dict(arrowstyle="-",lw=.55,color="#69727A"),bbox=dict(fc="white",ec="none",alpha=.72,pad=1))
    ax.set_title("FIG16  Frozen Cinch network: 167 loci, 223 edges",fontsize=17,weight="bold")
    ax.text(.01,.01,"Node size = degree · color = detected community · edge color = PP/PT/TP/TT · width = NMI",transform=ax.transAxes)
    save(fig,out/"FIG16_FINAL_NETWORK")
    # FIG17 community overview
    comm_rows=[]
    for i,c in enumerate(communities,1):
        sg=G.subgraph(c); top=sorted(sg.degree,key=lambda x:x[1],reverse=True)[:4]
        comm_rows.append({"community":i,"nodes":len(c),"edges":sg.number_of_edges(),"top_loci":";".join(n for n,_ in top),"top_annotations":";".join(str(G.nodes[n].get("annotation","")) for n,_ in top)})
    cdf=pd.DataFrame(comm_rows); cdf.to_csv(net/"NETWORK_COMMUNITIES.tsv",sep="\t",index=False)
    show=cdf.head(12)
    fig,ax=plt.subplots(figsize=(11,5.8)); x=np.arange(len(show)); ax.bar(x,show.nodes,color=[plt.cm.tab20((i-1)%20) for i in show.community]); ax.set_xticks(x,[f"C{i}" for i in show.community]); ax.set_ylabel("loci"); ax.set_title("FIG17  Connected-community architecture",fontsize=16,weight="bold")
    for xi,row in enumerate(show.itertuples()): ax.text(xi,row.nodes+1,str(row.top_loci).replace(";","\n"),ha="center",va="bottom",fontsize=7)
    ax.text(.01,.97,"Communities are descriptive graph partitions, not inferred pathways.",transform=ax.transAxes,va="top")
    fig.tight_layout(); save(fig,out/"FIG17_COMMUNITIES")
    return G,nd,ed


def copy_evidence(repo: Path, source: Path) -> None:
    diag=repo/"examples/spn534/diagnostics"; diag.mkdir(parents=True,exist_ok=True)
    mapping={
        source/"track_B_retrospective/figures/data/original_distance_diagnostic/DISTANCE_DEFINITIONS.md":diag/"LEGACY_DISTANCE_DEFINITIONS.md",
        source/"raw_type_channels/results/TABLE_SPNA_HIGH_PAIR_COINFINDER_SUMMARY.tsv":diag/"TABLE_COINFINDER_CONCORDANCE_SUMMARY.tsv",
        source/"raw_type_channels/results/TABLE_SPNA_HIGH_PAIR_ANNOTATION.tsv":diag/"TABLE_COINFINDER_PAIR_ANNOTATION.tsv",
        source/"track_B_retrospective/results/TABLE_B14_VATPASE_LINKAGE_CLASSIFICATION.tsv":diag/"TABLE_VATPASE_LINKAGE_CONTROL.tsv",
        source/"raw_type_channels/figures/data/SPN534_PP_VS_TT_TOP20_OVERLAY.tsv":repo/"examples/spn534/plot_data/PP_VS_TT_TOP20_OVERLAY.tsv",
        source/"final_neff_filtered/phandango/NEFF_FILTERED_PP_TT.tree.newick":repo/"examples/spn534/tree_state/FROZEN_PP_TT.tree.newick",
        source/"final_neff_filtered/phandango/NEFF_FILTERED_PP_TT.states.tsv.gz":repo/"examples/spn534/tree_state/FROZEN_PP_TT.states.tsv.gz",
        source/"final_neff_filtered/phandango/NEFF_FILTERED_PP_TT.columns.tsv":repo/"examples/spn534/tree_state/FROZEN_PP_TT.columns.tsv",
        source/"final_neff_filtered/phandango/NEFF_FILTERED_PP_TT.edges.tsv":repo/"examples/spn534/tree_state/FROZEN_PP_TT.edges.tsv",
    }
    for a,b in mapping.items(): b.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(a,b)
    pp_tt=pd.read_parquet(source/"raw_type_channels/figures/data/SPN534_PP_VS_TT_PLOTTING_DATA.parquet")
    pp_tt.to_csv(repo/"examples/spn534/plot_data/PP_VS_TT_PLOTTING_DATA.tsv.gz",sep="\t",index=False,compression="gzip")


def interactive_html(repo: Path, G: nx.Graph) -> None:
    net=repo/"examples/spn534/network"
    elements=[]
    for n,a in G.nodes(data=True):
        elements.append({"data":{"id":n,"label":n,"annotation":str(a.get("annotation","")),"degree":a["degree"],"community":a["community"]},"position":{"x":420*a["x"],"y":420*a["y"]}})
    for i,(u,v,a) in enumerate(G.edges(data=True)):
        data={"id":f"e{i}","source":u,"target":v,**a}; elements.append({"data":data})
    tpl='''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Cinch frozen SPN534 network</title><script src="https://unpkg.com/cytoscape@3.30.4/dist/cytoscape.min.js"></script><style>body{margin:0;font:14px system-ui;color:#17212b;background:#f6f8fa}header{padding:12px 18px;background:#16324f;color:white}#wrap{display:grid;grid-template-columns:300px 1fr 310px;height:calc(100vh - 55px)}aside{padding:14px;overflow:auto;background:white;border-right:1px solid #ddd}#detail{border-left:1px solid #ddd;border-right:0}#cy{min-height:600px}.control{margin:12px 0}.control label{display:block;font-weight:600}.chips label{display:inline-block;margin:3px;padding:4px 7px;background:#edf2f7;border-radius:12px;font-weight:500}button,input{max-width:100%}button{padding:6px 9px;margin:3px}.value{float:right;font-variant-numeric:tabular-nums}pre{white-space:pre-wrap;word-break:break-word;font-size:12px}</style></head><body><header><b>Cinch SPN534 frozen interaction network</b> — 167 loci / 223 edges</header><div id="wrap"><aside><div class="control"><label>Search locus/product</label><input id="search" placeholder="gyrB, transporter…"><button onclick="findNode()">Find</button><button onclick="resetView()">Reset</button></div><div class="control chips"><b>Channels</b><br><label><input type="checkbox" class="ch" value="PP" checked> PP</label><label><input type="checkbox" class="ch" value="PT" checked> PT</label><label><input type="checkbox" class="ch" value="TP" checked> TP</label><label><input type="checkbox" class="ch" value="TT" checked> TT</label><button onclick="crossOnly()">P↔T only</button></div><div id="sliders"></div><div class="control"><b>Edge color</b><select id="colorMode"><option value="channel">channel</option><option value="Neff">Neff</option><option value="HC_supported">HC support</option></select></div><p>Zoom, pan and drag directly. Click a node or edge for frozen metadata.</p><button onclick="spotlight('gyrB','aguA')">Highlight gyrB–aguA</button><button onclick="spotPP()">Four frozen PP edges</button></aside><main id="cy"></main><aside id="detail"><h3>Selection details</h3><pre id="info">Click an element.</pre></aside></div><script>const elements=__ELEMENTS__;const colors={PP:'#2676B8',PT:'#E69F00',TP:'#8B5CF6',TT:'#D84A4A'};const cy=cytoscape({container:document.getElementById('cy'),elements,layout:{name:'preset',fit:true,padding:30},style:[{selector:'node',style:{'background-color':'mapData(community,1,20,#2C73B9,#E07A5F)','width':'mapData(degree,1,62,10,44)','height':'mapData(degree,1,62,10,44)','label':'data(label)','font-size':7,'text-opacity':0.72}},{selector:'edge',style:{'line-color':e=>colors[e.data('channel')],'width':'mapData(NMI,0,1,0.5,5)','opacity':0.45,'curve-style':'bezier'}},{selector:'.dim',style:{opacity:.05}},{selector:'.focus',style:{'border-width':4,'border-color':'#fbbf24',opacity:1,'z-index':9999}},{selector:'edge.focus',style:{'line-color':'#ef4444','width':7,opacity:1}}]});
const defs=[['MIraw',0,2.3,0.01],['NMI',0,1,0.01],['Neff',0,10,0.1],['order_distance',0,1200,10],['HC_supported',0,70,1]];const s=document.getElementById('sliders');defs.forEach(([k,min,max,step])=>{s.insertAdjacentHTML('beforeend',`<div class="control"><label>${k} ≥ <span class="value" id="v_${k}">${min}</span></label><input id="s_${k}" type="range" min="${min}" max="${max}" step="${step}" value="${min}"></div>`);document.getElementById('s_'+k).oninput=applyFilters});document.querySelectorAll('.ch').forEach(x=>x.onchange=applyFilters);function applyFilters(){const allowed=new Set([...document.querySelectorAll('.ch:checked')].map(x=>x.value));cy.edges().forEach(e=>{let ok=allowed.has(e.data('channel'));defs.forEach(([k])=>{const v=+document.getElementById('s_'+k).value;document.getElementById('v_'+k).textContent=v;ok=ok&&(+e.data(k)>=v)});e.style('display',ok?'element':'none')});cy.nodes().forEach(n=>n.style('display',n.connectedEdges(':visible').length?'element':'none'))}function resetView(){cy.elements().removeClass('dim focus');document.querySelectorAll('.ch').forEach(x=>x.checked=true);defs.forEach(([k,min])=>document.getElementById('s_'+k).value=min);applyFilters();cy.fit(cy.elements(':visible'),30)}function findNode(){const q=document.getElementById('search').value.toLowerCase();const hit=cy.nodes().filter(n=>(n.id()+' '+n.data('annotation')).toLowerCase().includes(q));cy.elements().addClass('dim');hit.removeClass('dim').addClass('focus');hit.connectedEdges().removeClass('dim');if(hit.length)cy.animate({center:{eles:hit},zoom:2},{duration:400})}function spotlight(a,b){resetView();const ns=cy.getElementById(a).union(cy.getElementById(b));cy.elements().addClass('dim');ns.removeClass('dim').addClass('focus');ns.edgesWith(ns).removeClass('dim').addClass('focus');cy.fit(ns,140)}function spotPP(){resetView();const ee=cy.edges().filter(e=>e.data('channel')==='PP');cy.elements().addClass('dim');ee.removeClass('dim').addClass('focus');ee.connectedNodes().removeClass('dim').addClass('focus');cy.fit(ee.connectedNodes(),100)}function crossOnly(){document.querySelectorAll('.ch').forEach(x=>x.checked=['PT','TP'].includes(x.value));applyFilters()}cy.on('tap','node,edge',e=>{document.getElementById('info').textContent=JSON.stringify(e.target.data(),null,2)});document.getElementById('colorMode').onchange=e=>{const m=e.target.value;cy.edges().style('line-color',x=>m==='channel'?colors[x.data('channel')]:m==='Neff'?`hsl(${Math.min(220,20*x.data('Neff'))},70%,48%)`:`hsl(${Math.min(220,3*x.data('HC_supported'))},65%,45%)`)};</script></body></html>'''
    (net/"CINCH_INTERACTIVE.html").write_text(tpl.replace("__ELEMENTS__",json.dumps(elements,separators=(",",":"))),encoding="utf-8")


def sync_site_assets(repo: Path) -> None:
    site=repo/"site"; assets=site/"assets"; assets.mkdir(parents=True,exist_ok=True)
    figures=repo/"docs/assets/figures"
    for p in figures.glob("*.png"):
        shutil.copy2(p,assets/p.name)
    shutil.copy2(repo/"examples/spn534/network/CINCH_INTERACTIVE.html",site/"network.html")


def write_presentation_manifest(repo: Path) -> None:
    roots=[repo/"docs/assets/figures",repo/"examples/spn534/diagnostics",repo/"examples/spn534/network",repo/"examples/spn534/tree_state"]
    files=[p for r in roots for p in r.rglob("*") if p.is_file()]
    frozen=repo/"examples/spn534/results/FINAL_CANDIDATES.tsv"; files.append(frozen)
    rows=[]
    for p in sorted(set(files)):
        rows.append({"path":p.relative_to(repo).as_posix(),"bytes":p.stat().st_size,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()})
    pd.DataFrame(rows).to_csv(repo/"examples/spn534/results/PRESENTATION_MANIFEST.tsv",sep="\t",index=False)


def main() -> None:
    ap=argparse.ArgumentParser(); ap.add_argument("--repo",type=Path,default=Path(__file__).resolve().parents[1]); ap.add_argument("--source",type=Path,default=Path(r"D:\Codex_tools\Cinch\Cinch_full_handoff_20260825\external_validation\CINCH_EXT_SPN534")); args=ap.parse_args()
    repo=args.repo.resolve(); source=args.source.resolve(); out=repo/"docs/assets/figures"; out.mkdir(parents=True,exist_ok=True); style()
    frozen=repo/"examples/spn534/results/FINAL_CANDIDATES.tsv"; before=__import__("hashlib").sha256(frozen.read_bytes()).hexdigest()
    copy_evidence(repo,source); conceptual_figures(out); distance_figure(source,repo,out); data_figures(repo,source,out); G,_,_=build_network(repo,out); interactive_html(repo,G); sync_site_assets(repo); write_presentation_manifest(repo)
    after=__import__("hashlib").sha256(frozen.read_bytes()).hexdigest()
    assert before==after, "Frozen candidate table changed"
    print(json.dumps({"frozen_sha256":after,"edges":G.number_of_edges(),"nodes":G.number_of_nodes(),"figures":17},indent=2))


if __name__=="__main__": main()
