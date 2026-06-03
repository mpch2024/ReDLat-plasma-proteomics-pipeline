import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# =====================================================
# PATHS
# =====================================================

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RESULTS_DIR=rf"{PROJECT_DIR}\Graficos\DeLong"
os.makedirs(RESULTS_DIR,exist_ok=True)

DELONG_FILE=rf"{PROJECT_DIR}\Results_with_ptau\delong_results.csv"

# =====================================================
# LOAD
# =====================================================

results_delong=pd.read_csv(DELONG_FILE)

# =====================================================
# FILTER
# =====================================================

label_map={
"Soma":"Proteins",
"pTau":"p-tau217",
"Combo":"Proteins + p-tau217"
}

scenarios=["Proteins","p-tau217","Proteins + p-tau217"]

heatmap=pd.DataFrame(np.nan,index=scenarios,columns=scenarios)

for _,row in results_delong.iterrows():

    s1,s2=row[
    "Comparison"
    ].split(
    "_vs_"
    )

    s1=label_map[
    s1
    ]

    s2=label_map[
    s2
    ]

    heatmap.loc[
    s1,
    s2
    ]=row[
    "P"
    ]

    heatmap.loc[
    s2,
    s1
    ]=row[
    "P"
    ]

# =====================================================
# LABELS
# =====================================================

# mostrar solo significativos

mask=(

heatmap.isna()
|
(heatmap>0.05)
)
annot=heatmap.round(3).astype(str)
annot=annot.where(~mask,"")

# =====================================================
# FIGURE
# =====================================================

fig_w=6/2.54
fig_h=4/2.54
fig,ax=plt.subplots(

figsize=(fig_w,fig_h))

ax.set_position([
0.12,
0.15,
0.62,
0.72
])

hm=sns.heatmap(
heatmap,
cmap="viridis",
annot=annot,
fmt="",
mask=mask,
vmin=0,
vmax=0.005,
linewidths=0.5,
linecolor="white",
square=False,
annot_kws={"fontsize":4},
cbar_kws={
"label":"p-value",
"ticks":[0,0.001,0.002,0.003,0.004,0.005]},
ax=ax
)

# =====================================================
# COLORBAR
# =====================================================

cbar=hm.collections[0].colorbar

cbar.ax.tick_params(labelsize=5)

cbar.ax.yaxis.label.set_size(5)

# =====================================================
# STYLE
# =====================================================

ax.set_title(

"DeLong Test p-values",
fontsize=5,
pad=10
)

ax.set_xticklabels(
ax.get_xticklabels(),
rotation=45,
ha="right",
fontsize=5
)

ax.set_yticklabels(
ax.get_yticklabels(),
fontsize=5
)

# =====================================================
# SAVE
# =====================================================

plt.savefig(rf"{RESULTS_DIR}\DeLong_heatmap.pdf",
dpi=600,pad_inches=0)

plt.show()