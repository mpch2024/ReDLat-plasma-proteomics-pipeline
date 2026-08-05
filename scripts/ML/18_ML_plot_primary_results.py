###############################################################################
# ReDLat plasma proteomics — machine-learning workflow
# 18. Generate primary classification panels
# Requires: validated private analysis outputs
# Produces: editable panel files for Figure 4
# Data policy: participant-level inputs and predictions remain local.
###############################################################################

from pathlib import Path
import sys

_REPO_HINT = Path(__file__).resolve()
for _candidate in [_REPO_HINT.parent, *_REPO_HINT.parents]:
    if (_candidate / ".redlat-root").exists() or (_candidate / ".here").exists():
        PROJECT_ROOT = _candidate
        break
else:
    raise FileNotFoundError("Repository root not found. Set REDLAT_PROJECT_ROOT.")
sys.path.insert(0, str(PROJECT_ROOT / "python"))
from redlat_ml.config import load_config, require_files
CONFIG = load_config(__file__)

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from pathlib import Path
from matplotlib.patches import FancyBboxPatch

# =====================================================
# PATHS
# =====================================================

PROJECT_DIR = CONFIG.project_root
RESULTS_DIR = CONFIG.publication_root / "figures" / "figure_4" / "primary"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
ROC_FILE = CONFIG.private_root / "primary" / "02_nested_cv" / "roc_curves.csv"
PERM_FILE = CONFIG.private_root / "primary" / "02_nested_cv" / "Final_Panels" / "permutation_importance_freq08.csv"

# =====================================================
# LOAD
# =====================================================

roc = pd.read_csv(ROC_FILE)
perm = pd.read_csv(PERM_FILE)

# ---------------------------
# ROC curve chart
# ---------------------------

def plot_roc():

    data_size_cm=4.5

    left_cm=1.1
    bottom_cm=0.9
    right_cm=0.3
    top_cm=0.5

    fig_w=(left_cm+data_size_cm+right_cm)/2.54
    fig_h=(bottom_cm+data_size_cm+top_cm)/2.54

    fig,ax=plt.subplots(figsize=(fig_w,fig_h))

    ax.set_position([
    left_cm/(left_cm+data_size_cm+right_cm),
    bottom_cm/(bottom_cm+data_size_cm+top_cm),
    data_size_cm/(left_cm+data_size_cm+right_cm),
    data_size_cm/(bottom_cm+data_size_cm+top_cm)])

    colors={
    "Proteins":"#0057D9",
    #"p-tau217":"#D94B00",
    #"Proteins + p-tau217":"#009E73"
    "Proteins + APOE": "#F4B400"
    }

    text=[]

    mean_fpr=np.linspace(0,1,100)

    for model in roc["Model"].unique():
        tmp_model=roc[
        roc["Model"]==
        model
        ]
        tprs=[]
        aucs=[]
        for fold in sorted(
        tmp_model["Fold"].unique()
        ):
            tmp=tmp_model[
            tmp_model["Fold"]==
            fold
            ]
            interp=np.interp(
            mean_fpr,
            tmp["FPR"],
            tmp["TPR"]
            )

            interp[0]=0
            tprs.append(interp)
            aucs.append(
            tmp[
            "AUC"
            ].iloc[0])

        mean_tpr=np.mean(tprs,axis=0)
        sd_tpr=np.std(tprs,axis=0)
        mean_auc=np.mean(aucs)
        sd_auc=np.std(aucs)
        mean_tpr[-1]=1
        ax.plot(mean_fpr,mean_tpr,lw=1.3,color=colors[model])
        ax.fill_between(mean_fpr,np.clip(mean_tpr-sd_tpr,0,1),
        np.clip(mean_tpr+sd_tpr,0,1),color=colors[model],alpha=0.15)
        text.append((model,colors[model],f"{mean_auc:.2f}""\u00B1"f"{sd_auc:.2f}"))

    ax.plot([0,1],[0,1],"k--",lw=1,alpha=0.5)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_xticks([0,0.5,1])
    ax.set_yticks([0.5,1])
    ax.tick_params(axis="both",labelsize=6,length=0,pad=3)
    ax.set_xlabel("1 − Specificity",fontsize=7)
    ax.set_ylabel("Sensitivity",fontsize=7)
    
    # ---------- legend box ----------

    x0 = 0.19
    y0 = 0.03

    box_w = 0.78
    line_h = 0.05

    box_h = 0.02 + len(text)*line_h
         
    rect = FancyBboxPatch(
    (x0, y0),
    box_w,
    box_h,
    boxstyle="round,pad=0.015",
    facecolor="#F8F8F8",
    edgecolor="#777777",
    #linewidth=0.45,
    transform=ax.transAxes,
    zorder=1
    )  
    
    ax.add_patch(rect)

    for i,(name,color,val) in enumerate(text):

        yy = y0 + box_h - (i+1)*line_h + 0.012

        ax.plot(
        [x0+0.03, x0+0.08],
        [yy, yy],
        color=color,
        lw=1.6,
        solid_capstyle="round",
        transform=ax.transAxes,
        zorder=3
        )

        ax.text(
        x0+0.10,
        yy,
        f"{name}: {val}",
        fontsize=5.3,
        color="#303030",
        ha="left",
        va="center",
        transform=ax.transAxes,
        zorder=3,
        clip_on=False
        )

# ---------- end legend ---------- 

    for s in ax.spines.values():
        s.set_linewidth(0.6)
    plt.savefig(RESULTS_DIR/"ROC_curve.pdf",dpi=600,pad_inches=0)

    plt.show()

# =====================================================
# PERMUTATION
# =====================================================

def plot_permutation():

    data_size_cm = 4.5
    left_cm = 2.0
    bottom_cm = 0.9
    right_cm = 0.3
    top_cm = 0.5

    fig_w = (left_cm + data_size_cm + right_cm) / 2.54
    fig_h = (bottom_cm + data_size_cm + top_cm) / 2.54
    fig,ax = plt.subplots(figsize=(fig_w,fig_h))

    ax.set_position([
        left_cm /(left_cm + data_size_cm + right_cm),
        bottom_cm /(bottom_cm + data_size_cm + top_cm),
        data_size_cm /(left_cm + data_size_cm + right_cm),
        data_size_cm /(bottom_cm + data_size_cm + top_cm)
    ])

    perm_plot = (perm.sort_values("Permutation",ascending=True))

    ax.barh(
    perm_plot["Proteins"],
    perm_plot["Permutation"],
    xerr=perm_plot["Permutation_sd"],
    color="#F4B400", #cambiar color
    linewidth=0.5
    )
    ax.margins(y=0.02)

    ax.set_xlabel("Feature Importance",fontsize=7)
    ax.set_ylabel("")
    ax.tick_params(
        axis="x",
        labelsize=7,
        pad=0,
        length=2,
        width=0.5
    )

    ax.tick_params(
        axis="y",
        labelsize=6,
        pad=1,
        length=2,
        width=0.5
    )

    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.set_title("Proteins",fontsize=7)
    plt.savefig(RESULTS_DIR/"Permutation_importance.pdf",dpi=600,pad_inches=0)
    plt.show()

# ===========================
# RUN
# ===========================

plot_roc()
plot_permutation()  

