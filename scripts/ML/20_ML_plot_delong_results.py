###############################################################################
# ReDLat plasma proteomics — machine-learning workflow
# 20. Generate model-comparison panel
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
import seaborn as sns
from pathlib import Path

# =====================================================
# PATHS
# =====================================================

PROJECT_DIR = CONFIG.project_root
RESULTS_DIR = CONFIG.publication_root / "figures" / "figure_4" / "delong"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
DELONG_FILE = CONFIG.private_root / "fixed_panel_ptau" / "delong_results.csv"

# LOAD
results_delong=pd.read_csv(DELONG_FILE)

# =====================================================
# FILTER
# =====================================================

label_map={
#"Soma":"Proteins",
"Panel":"Proteins",
"pTau":"p-tau217",
"Combo":"Combined"
}

scenarios=["Proteins","p-tau217","Combined"]

heatmap=pd.DataFrame(np.nan,index=scenarios,columns=scenarios)

for _,row in results_delong.iterrows():

    s1,s2=row["Comparison"].split("_vs_")
    s1=label_map[s1]
    s2=label_map[s2]
    heatmap.loc[s1,s2]=row["P"]
    heatmap.loc[s2,s1]=row["P"]

# LABELS (mostrar solo significativos)
mask=(heatmap.isna()|(heatmap>0.05))
annot=heatmap.round(3).astype(str)
annot=annot.where(~mask,"")

# =====================================================
# FIGURE
# =====================================================

fig_w=5/2.54
fig_h=4/2.54
fig,ax=plt.subplots(figsize=(fig_w,fig_h))

plt.subplots_adjust(
    left=0.25,
    bottom=0.25,
    right=0.90,
    top=0.90
)

hm=sns.heatmap(
heatmap,
cmap="viridis",
annot=annot,
fmt="",
mask=mask,
vmin=0,
vmax=0.05,
linewidths=0.5,
linecolor="white",
square=False,
annot_kws={"fontsize":4},
cbar_kws={
"label":"p-value",
"ticks":[0,0.01,0.02,0.03,0.04,0.05]},
ax=ax
)

# COLORBAR
cbar=hm.collections[0].colorbar
cbar.ax.tick_params(labelsize=5)
cbar.ax.yaxis.label.set_size(5)

# STYLE
ax.set_title("DeLong Test p-values",fontsize=5,pad=5)
ax.set_xticklabels(ax.get_xticklabels(),rotation=45,ha="right",fontsize=5)
ax.set_yticklabels(ax.get_yticklabels(),fontsize=5)

# SAVE
plt.savefig(RESULTS_DIR/"DeLong_heatmap_panel_ptau.pdf",dpi=600,pad_inches=0)

plt.show()