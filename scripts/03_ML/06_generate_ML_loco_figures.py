import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from matplotlib.patches import FancyBboxPatch

# =====================================================
# PATHS
# =====================================================

PROJECT_DIR = r"C:\Clasificacion"
RESULTS_DIR = rf"{PROJECT_DIR}\Graficos\LOCO"
ROC_FILE = rf"{PROJECT_DIR}\Results_LOCO\roc_curves.csv"
ROC_FILE_POOLED = rf"{PROJECT_DIR}\Results\roc_curves.csv"
PERM_FILE = rf"{PROJECT_DIR}\Results_LOCO\permutation_importance.csv"
PERM_FILE_COUNTRY = rf"{PROJECT_DIR}\Results_LOCO\permutation_importance_country.csv"

# =====================================================
# LOAD
# =====================================================

roc = pd.read_csv(ROC_FILE)
roc_pooled = pd.read_csv(ROC_FILE_POOLED)
perm = pd.read_csv(PERM_FILE)
perm_country = pd.read_csv(PERM_FILE_COUNTRY)
# ---------------------------
# ROC curve chart
# ---------------------------

def plot_roc():

    data_size_cm=5
    fig,ax=plt.subplots(figsize=(2.3,2.3))
    
    colors={
    "Argentina":"#795548",
    "Chile":"#8B0000",
    "Colombia":"#56B4E9",
    "Mexico":"#8E7CC3",
    "Peru":"#4D4D4D"
    }

    text=[]

    mean_fpr=np.linspace(0,1,100)

    for country in roc["Country"].unique():

        tmp=roc[roc["Country"]==country]
        interp=np.interp(mean_fpr,tmp["FPR"],tmp["TPR"])
        interp[0]=0
        auc_country=tmp[
        "AUC"].iloc[0]

        ax.plot(mean_fpr,interp,lw=1.4,color=colors[country])
        text.append((country,colors[country],f"{auc_country:.2f}"))

    ax.plot([0,1],[0,1],"k--",lw=1,alpha=0.5)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_xticks([0,0.5,1])
    ax.set_yticks([0,0.5,1])
    ax.tick_params(labelsize=6,pad=3,length=0)
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
        f"{name} (AUC={val})",
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
    plt.savefig(rf"{RESULTS_DIR}\01_ROC_country.pdf",dpi=600,pad_inches=0)

    plt.show()


def plot_roc_global():

    fig,ax=plt.subplots(figsize=(2.3,2.3))

    mean_fpr=np.linspace(0,1,100)
    models=[
    (roc,"LOCO","#D94B00"),
    (roc_pooled,"Pooled","#0057D9")]

    text=[]

    for data,label,color in models:
        tprs=[]
        aucs=[]

        for fold in sorted(
        data[
        "Fold"
        ].unique()
        ):
            tmp=data[data["Fold"]==fold]
            interp=np.interp(
            mean_fpr,tmp["FPR"],
            tmp["TPR"])
            interp[0]=0
            tprs.append(interp)
            aucs.append(
            tmp["AUC"].iloc[0])

        mean_tpr=np.mean(tprs,axis=0)
        sd=np.std(tprs,axis=0)
        mean_auc=np.mean(aucs)
        sd_auc=np.std(aucs)

        ax.plot(
        mean_fpr,
        mean_tpr,
        lw=1.5,
        color=color
        )

        ax.fill_between(
        mean_fpr,
        np.clip(mean_tpr-sd,0,1),
        np.clip(mean_tpr+sd,0,1),
        color=color,
        alpha=0.15
        )

        text.append((label,color,
        f"{mean_auc:.2f}±{sd_auc:.2f}"))

    ax.plot([0,1],[0,1],"k--",alpha=0.5)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_xticks([0,0.5,1])
    ax.set_yticks([0,0.5,1])
    ax.tick_params(labelsize=6,length=0,pad=3)
    ax.set_xlabel("1 − Specificity",fontsize=7)
    ax.set_ylabel("Sensitivity",fontsize=7)
    
    box=FancyBboxPatch(
    (0.38,0.03),
    0.58,
    0.12,
    boxstyle="round,pad=0.015",
    facecolor="#F8F8F8",
    edgecolor="#777777",
    transform=ax.transAxes
    )

    ax.add_patch(box)

    for i,(name,color,val) in enumerate(text):

        yy=0.11-i*0.05

        ax.plot(
        [0.40,0.46],
        [yy,yy],
        color=color,
        lw=1.5,
        transform=ax.transAxes
        )

        ax.text(
        0.50,
        yy,
        f"{name}: {val}",
        fontsize=5.5,
        transform=ax.transAxes,
        va="center"
        )

    plt.savefig(rf"{RESULTS_DIR}\ROC_global.pdf",dpi=600,pad_inches=0)
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
    perm_plot["Protein"],
    perm_plot["Permutation"],
    xerr=perm_plot["Permutation_sd"],
    color="#F4B400",
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
    plt.savefig(rf"{RESULTS_DIR}\Permutation_importance.pdf",dpi=600,pad_inches=0)
    plt.show()

# =====================================================
# PERMUTATION COUNTRY
# =====================================================

country_colors={
"Argentina":"#795548",
"Chile":"#8B0000",
"Colombia":"#56B4E9",
"Mexico":"#8E7CC3",
"Peru":"#4D4D4D"
}


def plot_permutation_country():

    countries=list(
    perm_country[
    "Country"
    ].unique()
    )

    for country in countries:

        tmp=(
        perm_country[
        perm_country[
        "Country"
        ]==
        country
        ]
        .sort_values(
        "Permutation_mean",
        ascending=True))

        data_size_cm=4.5

        fig,ax=plt.subplots(figsize=(2.5,2.3))

        ax.barh(
        tmp["Protein"],
        tmp["Permutation_mean"],
        xerr=tmp["Permutation_sd"],
        color=country_colors[country],
        linewidth=0.5
        )

        ax.set_xlabel("Feature Importance",fontsize=7)
        ax.set_ylabel("")
        ax.tick_params(
        axis="x",
        labelsize=6,
        length=0,
        pad=2
        )

        ax.tick_params(axis="y",labelsize=6,length=0,pad=1)
        ax.set_title(country,fontsize=7)

        for s in ax.spines.values():
            s.set_linewidth(0.5)

        plt.savefig(rf"{RESULTS_DIR}\Permutation_{country}.pdf",
        dpi=600,pad_inches=0)
        plt.close()

# ===========================
# RUN
# ===========================

plot_roc()
plot_permutation()   
plot_roc_global()
plot_permutation_country()

#0057D9: azul
#D94B00: rojo anaranjado
#009E73: verde
#F4B400: amarillo