###############################################################################
# ReDLat plasma proteomics — machine-learning workflow
# 12. Run country-held-out classification
# Requires: Scripts 10–11 outputs
# Produces: country-specific performance and private OOF predictions
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
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    auc,
    roc_auc_score,
    roc_curve,
    accuracy_score,
    recall_score,
    precision_score,
    f1_score,
    confusion_matrix
)

from sklearn.linear_model import LogisticRegressionCV
from sklearn.svm import SVC
from sklearn.inspection import permutation_importance
from skopt import BayesSearchCV
from skopt.space import Real, Categorical
import warnings
warnings.simplefilter("ignore", FutureWarning)

# -----------------------------------------------------
# REPRODUCIBILITY
# -----------------------------------------------------

RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

# -----------------------------------------------------
# PATHS
# -----------------------------------------------------
PROJECT_DIR = CONFIG.project_root
DATA_DIR = CONFIG.data_dir
RESULTS_DIR = CONFIG.private_root / "loco" / "02_nested_cv"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FOLDS_DIR = CONFIG.private_root / "loco" / "00_folds"
DEP_DIR = CONFIG.private_root / "loco" / "01_dep_folds"
DATA_FILE = CONFIG.master_file

# -----------------------------------------------------
# PARAMETERS
# -----------------------------------------------------

TARGET = "SampleGroup"
ID_COL = "SampleId"
OUTER_FOLDS = 5
INNER_FOLDS = 5
MAX_ITER = 50000
TOL = 1e-3

# -----------------------------------------------------
# LOAD/Identificadores
# -----------------------------------------------------

df = pd.read_csv(DATA_FILE)

META_COLS = [
    "SampleGroup",
    "Site",
    "Country",
    "Sex",
    "Age",
    "Education",
    "ApoE",
    "APOE_group",
    "APOE4_carrier",
    "cdr_global",
    "cdr_boxscore",
    "mmse_total",
    "udsfaq_total",
    "cog_benson",
    "cog_tmt_a",
    "cog_tmt_b",
    "cog_craft_verb_delayed",
    "NPI",
    "Mini.SEA",
    "T.ADLQ",
    "p.tau217",
    "p.tau181",
    "NfL",
    "ratio.AB42.40",
    "GFAP_1"
]

df[ID_COL] = df[ID_COL].astype(str)
df = df.set_index(ID_COL)
df = df[df[TARGET].isin(["CN","AD"])].copy()
protein_cols = [
    c
    for c in df.columns
    if c not in META_COLS
]

print("\nMeta:",len(META_COLS))
print("\nProteins:",len(protein_cols))

y = df[TARGET].map({"CN":0,"AD":1})
X = df[protein_cols]

# =====================================================
# COUNTRY
# =====================================================

site_to_country = {
'BN':'Argentina',
'BE':'Chile',
'SL':'Chile',
'MA':'Colombia',
'LO':'Colombia',
'AF':'Mexico',
'CU':'Peru'
}

# Reviewer IDs are pseudonymous; country comes from released metadata.
df["country"] = df["Country"]
countries = sorted(df["country"].dropna().astype(str).unique().tolist())

print("\nCountries")
print(countries)

# =====================================================
# SVM SEARCH SPACE
# =====================================================

svm_param_space = {
    "kernel": Categorical(["linear","rbf"]),
    "C": Real(0.1,100,prior="log-uniform"),
    "gamma": Real(1e-4,1,prior="log-uniform")
}

# =====================================================
# LOCO
# =====================================================

metrics=[]
selection_summary=[]
outer_panels=[]
best_params=[]
roc_store=[]
INNER_STABILITY = 1.0

pipeline_params = pd.DataFrame([{
    "OUTER_FOLDS": OUTER_FOLDS,
    "INNER_FOLDS": INNER_FOLDS,
    "INNER_STABILITY": INNER_STABILITY,
    "MAX_ITER": MAX_ITER,
    "TOL": TOL,
    "RANDOM_STATE": RANDOM_STATE,
    "Cs": "np.logspace(-4,0,40)",
    "l1_ratios": "0.85,0.90,0.95,1.0",
    "SVM_iter": 50
}])

# =====================================================
# LOOP INNER CV → STABILITY
# =====================================================

oof = pd.Series(np.nan,index=df.index)

for country in countries:

    train_ids = pd.read_csv(FOLDS_DIR / f"{country}_train_ids.csv")
    test_ids = pd.read_csv(FOLDS_DIR / f"{country}_test_ids.csv")
    train_ids = train_ids["SampleId"].astype(str)
    test_ids  = test_ids["SampleId"].astype(str)
    X_outer_train = X.loc[train_ids].copy()
    X_outer_test  = X.loc[test_ids].copy()
    y_outer_train = y.loc[train_ids]
    y_outer_test  = y.loc[test_ids]

    dep_file = (DEP_DIR/country/"candidate_gene_symbols.csv")
    dep_genes = pd.read_csv(dep_file)
    dep_genes = (
        dep_genes["EntrezGeneSymbol"]
        .dropna()
        .astype(str)
        .tolist()
    )

    print(f"\nCountry {country}")
    print(f"DEP genes file: {len(dep_genes)}")
    missing = sorted(set(dep_genes) - set(X.columns))
    print(f"Missing genes: {len(missing)}")
    if len(missing) > 0:
        print(missing)

    keep = np.array([
        g
        for g in dep_genes
        if g in X.columns
    ])

    print(f"Genes retained: {len(keep)}")

    X_outer_train = X_outer_train[keep].copy()
    X_outer_test  = X_outer_test[keep].copy()

    print(f"Country {country} | DEP proteins: {len(keep)}")

    inner_selected=[]
    inner_cv=StratifiedKFold(
        n_splits=INNER_FOLDS,
        shuffle=True,
        random_state=RANDOM_STATE
    )

    for itr,ival in inner_cv.split(X_outer_train,y_outer_train):

        Xin=X_outer_train.iloc[itr]
        yin=y_outer_train.iloc[itr]
        scaler=StandardScaler()
        Xin=scaler.fit_transform(Xin)
        selector=LogisticRegressionCV(
            penalty="elasticnet",
            solver="saga",
            Cs=np.logspace(-4,0,40),
            l1_ratios=[0.85,0.90,0.95,1.0],
            cv=INNER_FOLDS,
            scoring="roc_auc",
            max_iter=MAX_ITER,
            tol=TOL,
            n_jobs=-1,
            random_state=RANDOM_STATE,
        )

        selector.fit(Xin,yin)
        coef=np.squeeze(selector.coef_)
        mask=np.abs(coef)>5e-3
        inner_selected.append(np.array(keep)[mask])

    # -----------------------------------------
    # CONSENSUS PANEL
    # -----------------------------------------

    freq={}
    for s in inner_selected:
        for p in s:
            freq[p]=freq.get(p,0)+1

    freq=pd.Series(freq)/INNER_FOLDS

    panel=freq[freq>=INNER_STABILITY].index.tolist()

    if len(panel)==0:
        panel=freq.sort_values(ascending=False).head(10).index.tolist()

    outer_panels.append(panel)
    pd.Series(panel).to_csv(RESULTS_DIR/f"panel_fold_{country}.csv",index=False)
    n_selected=len(panel)
    selection_summary.append({"Country":country,"N_features":n_selected})

    print("\nFrequency summary")
    print(freq.describe())
    print("\nTop stable proteins")
    print(freq.sort_values(ascending=False).head(10))
    print(f"\nConsensus panel: {n_selected}")

    # -----------------------------------------
    # TRAIN SOMASCAN SVM
    # -----------------------------------------

    X_train=X_outer_train[panel]
    X_test=X_outer_test[panel]

    scaler=StandardScaler()

    X_train=scaler.fit_transform(X_train)
    X_test=scaler.transform(X_test)

    search=BayesSearchCV(
        estimator=SVC(probability=False,random_state=RANDOM_STATE),
        search_spaces=svm_param_space,
        n_iter=50,
        cv=INNER_FOLDS,
        scoring="roc_auc",
        n_jobs=-1,
        random_state=RANDOM_STATE
    )

    search.fit(X_train,y_outer_train)
    y_prob=search.best_estimator_.decision_function(X_test)
    oof[test_ids]=(y_prob-y_prob.mean())/(y_prob.std()+1e-8)
    y_pred=(y_prob>=0).astype(int)
    tn,fp,fn,tp=confusion_matrix(y_outer_test,y_pred).ravel()
    auc=roc_auc_score(y_outer_test,y_prob)
    best_params.append({
        "Country": country,
        "kernel": search.best_params_["kernel"],
        "C": search.best_params_["C"],
        "gamma": search.best_params_["gamma"],
        "DEP_proteins": len(keep),
        "Consensus_panel": n_selected,
        "AUC": auc
    })

    fpr,tpr,thr=roc_curve(y_outer_test,y_prob)
    roc_store.append(pd.DataFrame({
    "Country":country,"Model":"Proteins","FPR":fpr,"TPR":tpr,"Threshold":thr,"AUC":auc}))

    metrics.append({
        "Country":country,
        "AUC":auc,
        "Accuracy":accuracy_score(y_outer_test,y_pred),
        "Sensitivity":recall_score(y_outer_test,y_pred),
        "Specificity":tn/(tn+fp),
        "Precision":precision_score(y_outer_test,y_pred,zero_division=np.nan),
        "F1":f1_score(y_outer_test,y_pred,zero_division=np.nan)
    })
    print(f"AUC={auc:.3f}")

# =====================================================
# SAVE
# =====================================================

metrics_df = pd.DataFrame(metrics)
metrics_df.to_csv(RESULTS_DIR/"metrics_nestedCV.csv",index=False)
pd.DataFrame(best_params).to_csv(RESULTS_DIR/"best_params_outer.csv",index=False)
pd.DataFrame(selection_summary).to_csv(RESULTS_DIR/"feature_counts.csv",index=False)
pipeline_params.to_csv(RESULTS_DIR/"pipeline_parameters.csv",index=False)

# ==========================================
# PANEL STABILITY
# ==========================================

panel_freq={}

for panel in outer_panels:
    for p in panel:
        panel_freq[p]=panel_freq.get(p,0)+1

panel_freq=(pd.Series(panel_freq)/len(outer_panels))
panel_freq=panel_freq.sort_values(ascending=False)
panel_freq.to_csv(RESULTS_DIR/"outer_panel_frequency.csv")

print("\n=======================")
print(metrics_df.mean(numeric_only=True))

metrics_df.mean(numeric_only=True).to_csv(RESULTS_DIR/"metrics_mean.csv")

# =====================================================
# OOF EXPORT
# =====================================================

oof_df=pd.DataFrame({"SampleId":df.index,"y_true":y,"oof":oof})
oof_df.to_csv(RESULTS_DIR/"oof_predictions.csv",index=False)

# =====================================================
# ROC
# =====================================================

roc_df=pd.concat(roc_store,ignore_index=True)
roc_df.to_csv(RESULTS_DIR/"roc_curves.csv",index=False)

# =====================================================
# SCORE SUMMARY
# =====================================================

score_summary=pd.DataFrame({
"Model":["SomaScan"],"Mean_CN":[oof[y==0].mean()],"Mean_AD":[oof[y==1].mean()]})
score_summary.to_csv(RESULTS_DIR/"score_summary.csv",index=False)

# =====================================================
# FINAL MODEL + PERMUTATION
# =====================================================

FINAL_PANEL_DIR = RESULTS_DIR / "Final_Panels"
FINAL_PANEL_DIR.mkdir(exist_ok=True)

final_panels = {
    "freq06": panel_freq[panel_freq >= 0.6].index.tolist(),
    "freq08": panel_freq[panel_freq >= 0.8].index.tolist(),
    "freq10": panel_freq[panel_freq >= 1.0].index.tolist()
}

for tag, panel_final in final_panels.items():

    if len(panel_final)==0:
        panel_final=(panel_freq.head(10).index.tolist())

    pd.Series(panel_final).to_csv(FINAL_PANEL_DIR/f"panel_final_{tag}.csv",index=False)
    best=pd.DataFrame(best_params)
    kernel=best["kernel"].mode()[0]
    C=best["C"].median()
    gamma=best["gamma"].median()
    X_final=X[panel_final].copy()
    scaler=StandardScaler()
    X_final=scaler.fit_transform(X_final)
    svm=SVC(
        kernel=kernel,
        C=C,
        gamma=gamma,
        probability=True,
        random_state=RANDOM_STATE
    )

    svm.fit(X_final,y)

    perm=permutation_importance(
        svm,
        X_final,
        y,
        n_repeats=100,
        scoring="roc_auc",
        random_state=RANDOM_STATE,
        n_jobs=-1
    )

    perm_summary=pd.DataFrame({
        "Proteins":panel_final,
        "Permutation":perm.importances_mean,
        "Permutation_sd":perm.importances_std
    })

    perm_summary=perm_summary.sort_values("Permutation",ascending=False)
    perm_summary.to_csv(FINAL_PANEL_DIR/f"permutation_importance_{tag}.csv",index=False)
