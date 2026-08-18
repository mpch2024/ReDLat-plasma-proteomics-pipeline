###############################################################################
# ReDLat plasma proteomics — machine-learning workflow
# 09. Evaluate APOE model extension
# Requires: Scripts 07–08 outputs
# Produces: nested model comparison and DeLong test
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
RESULTS_DIR = CONFIG.private_root / "apoe" / "02_nested_cv"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FOLDS_DIR = CONFIG.private_root / "apoe" / "00_folds"
DEP_DIR = CONFIG.private_root / "apoe" / "01_dep_folds"
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
# LOAD
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

df[ID_COL]=df[ID_COL].astype(str)
df=df.set_index(ID_COL)
df=df[df[TARGET].isin(["CN","AD"])].copy()

protein_cols=[
    c
    for c in df.columns
    if c not in META_COLS
]

print("\nMeta:",len(META_COLS))
print("\nProteins:",len(protein_cols))

y = df[TARGET].map({"CN":0,"AD":1})
X = df[protein_cols]

# =====================================================
# COMMON COHORT
# proteins + apoe only
# =====================================================

N0=len(df)
COMMON_COLS=protein_cols+["APOE4_carrier"]
df=df.dropna(subset=COMMON_COLS).copy()

print(df.index.name)
print(df.columns.tolist()[:10])

N1=len(df)

print("\n=======================")
print(f"Original N={N0}")
print(f"Final N={N1}")
print(f"Removed={N0-N1}")
print("=======================\n")

pd.DataFrame({
    "Original_N":[N0],
    "Final_N":[N1],
    "Removed":[N0-N1]
}).to_csv(RESULTS_DIR/"selected_subjects.csv",index=False)

sample_ids=df.index.astype(str)

y=df[TARGET].map({"CN":0,"AD":1})

X_soma=df[protein_cols]
X_apoe=df[["APOE4_carrier"]]
X = X_soma.copy()

print("\nSoma:",X_soma.shape)
print("\napoe:",X_apoe.shape)

cohort_export = df[[TARGET, "APOE4_carrier"]].reset_index()
cohort_export.to_csv(RESULTS_DIR/"cohort_common.csv",index=False)

# =====================================================
# SVM SEARCH SPACE
# =====================================================

svm_param_space = {
    "kernel": Categorical(["linear","rbf"]),
    "C": Real(0.1,100,prior="log-uniform"),
    "gamma": Real(1e-4,1,prior="log-uniform")
}

# =====================================================
# OUTER CV
# =====================================================

metrics_soma=[]
metrics_combo=[]
outer_panels=[]
roc_store=[]

oof_soma=np.zeros(len(df))
oof_combo=np.zeros(len(df))

best_params_soma=[]
best_params_combo=[]

selection_summary_soma=[]

INNER_STABILITY = 1.0

# =====================================================
# LOOP INNER CV → STABILITY
# =====================================================
print(df.index[:5])
print(X_soma.index[:5])
print(y.index[:5])

oof_soma = pd.Series(np.nan,index=df.index.astype(str))
oof_combo = pd.Series(np.nan,index=df.index.astype(str))

for fold in range(1,6):

    print(f"\nFold {fold}")

    train_ids = pd.read_csv(FOLDS_DIR / f"fold_{fold}_train_ids.csv")
    test_ids = pd.read_csv(FOLDS_DIR / f"fold_{fold}_test_ids.csv")
    train_ids = train_ids["SampleId"].astype(str)
    test_ids  = test_ids["SampleId"].astype(str)
    X_outer_train = X.loc[train_ids].copy()
    X_outer_test  = X.loc[test_ids].copy()
    y_outer_train = y.loc[train_ids]
    y_outer_test  = y.loc[test_ids]

    dep_file = (DEP_DIR/f"fold_{fold}"/"candidate_gene_symbols.csv")
    dep_genes = pd.read_csv(dep_file)
    dep_genes = (
        dep_genes["EntrezGeneSymbol"]
        .dropna()
        .astype(str)
        .tolist()
    )

    print(f"\nFold {fold}")
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

    print(f"Fold {fold} | DEP proteins: {len(keep)}")

    inner_selected=[]

    print(X_outer_train.shape)
    print(X_outer_test.shape)

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
    freq.sort_values(ascending=False).to_csv(RESULTS_DIR/f"frequency_fold_{fold}.csv")
    panel=freq[freq>=INNER_STABILITY].index.tolist()

    if len(panel)==0:
        panel=freq.sort_values(ascending=False).head(10).index.tolist()

    outer_panels.append(panel)
    pd.Series(panel).to_csv(RESULTS_DIR/f"panel_fold_{fold}.csv",index=False)
    n_selected=len(panel)
    selection_summary_soma.append({"Fold":fold,"N_features":n_selected})

    print(f"\nConsensus panel: {n_selected}")

    # =================================================
    # TRAIN SOMASCAN SVM
    # =================================================
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
    oof_soma.loc[test_ids] = (y_prob-y_prob.mean())/(y_prob.std()+1e-8)
    y_pred=(y_prob>=0).astype(int)
    tn,fp,fn,tp=confusion_matrix(y_outer_test,y_pred).ravel()
    auc = roc_auc_score(y_outer_test,y_prob)
    best_params_soma.append({
        "Fold": fold,
        "kernel": search.best_params_["kernel"],
        "C": search.best_params_["C"],
        "gamma": search.best_params_["gamma"],
        "DEP_proteins": len(keep),
        "Consensus_panel": n_selected,
        "AUC": auc
    })

    fpr,tpr,thr=roc_curve(y_outer_test,y_prob)
    roc_store.append(pd.DataFrame({
    "Fold":fold,"Model":"Proteins","FPR":fpr,"TPR":tpr,"Threshold":thr,"AUC":auc}))

    metrics_soma.append({
        "Fold":fold,
        "AUC":auc,
        "Accuracy":accuracy_score(y_outer_test,y_pred),
        "Sensitivity":recall_score(y_outer_test,y_pred),
        "Specificity":tn/(tn+fp),
        "Precision":precision_score(y_outer_test,y_pred,zero_division=np.nan),
        "F1":f1_score(y_outer_test,y_pred,zero_division=np.nan)
    })
    print(f"AUC={auc:.3f}")

    # =================================================
    # SOMASCAN + PTAU217
    # =================================================

    X_train_combo=pd.concat([
    X_outer_train[panel].reset_index(drop=True),
    X_apoe.loc[train_ids].reset_index(drop=True)
    ],axis=1)

    X_test_combo=pd.concat([
    X_outer_test[panel].reset_index(drop=True),
    X_apoe.loc[test_ids].reset_index(drop=True)
    ],axis=1)

    scaler=StandardScaler()

    X_train_combo=scaler.fit_transform(X_train_combo)
    X_test_combo=scaler.transform(X_test_combo)

    search=BayesSearchCV(
        estimator=SVC(probability=False,random_state=RANDOM_STATE),
        search_spaces=svm_param_space,
        n_iter=50,
        cv=INNER_FOLDS,
        scoring="roc_auc",
        n_jobs=-1,
        random_state=RANDOM_STATE
    )

    search.fit(X_train_combo,y_outer_train)
    y_prob=search.best_estimator_.decision_function(X_test_combo)
    oof_combo.loc[test_ids]=(y_prob-y_prob.mean())/(y_prob.std()+1e-8)
    y_pred=(y_prob>=0).astype(int)
    tn,fp,fn,tp=confusion_matrix(y_outer_test,y_pred).ravel()
    auc=roc_auc_score(y_outer_test,y_prob)
    best_params_combo.append({
        "Fold": fold,
        "kernel": search.best_params_["kernel"],
        "C": search.best_params_["C"],
        "gamma": search.best_params_["gamma"],
        "DEP_proteins": len(keep),
        "Consensus_panel": n_selected,
        "AUC": auc
    })

    fpr,tpr,thr=roc_curve(y_outer_test,y_prob)
    roc_store.append(pd.DataFrame({
    "Fold":fold,"Model":"Proteins + APOE","FPR":fpr,"TPR":tpr,"Threshold":thr,"AUC":auc}))

    metrics_combo.append({
    "Fold":fold,
    "AUC":auc,
    "Accuracy":accuracy_score(y_outer_test,y_pred),
    "Sensitivity":recall_score(y_outer_test,y_pred),
    "Specificity":tn/(tn+fp),
    "Precision":precision_score(y_outer_test,y_pred),
    "F1":f1_score(y_outer_test,y_pred)
    })

    print(f"Combo AUC={metrics_combo[-1]['AUC']:.3f}")

# =====================================================
# SAVE METRICS
# =====================================================

metrics_soma=pd.DataFrame(metrics_soma)
metrics_combo=pd.DataFrame(metrics_combo)

metrics_soma.to_csv(RESULTS_DIR/"metrics_soma.csv",index=False)
metrics_combo.to_csv(RESULTS_DIR/"metrics_combo.csv",index=False)

pd.DataFrame(best_params_soma).to_csv(RESULTS_DIR/"best_params_outer_soma.csv",index=False)
pd.DataFrame(best_params_combo).to_csv(RESULTS_DIR/"best_params_outer_combo.csv",index=False)

pd.DataFrame(selection_summary_soma).to_csv(RESULTS_DIR/"feature_counts_soma.csv",index=False)


# =====================================================
# MODEL COMPARISON SUMMARY
# =====================================================

summary = pd.DataFrame({
    "Model": ["SomaScan", "Combo"],
    "Mean_AUC": [
        metrics_soma["AUC"].mean(),
        metrics_combo["AUC"].mean()
    ],
    "SD_AUC": [
        metrics_soma["AUC"].std(),
        metrics_combo["AUC"].std()
    ],
    "Mean_Accuracy": [
        metrics_soma["Accuracy"].mean(),
        metrics_combo["Accuracy"].mean()
    ],
    "Mean_Sensitivity": [
        metrics_soma["Sensitivity"].mean(),
        metrics_combo["Sensitivity"].mean()
    ],
    "Mean_Specificity": [
        metrics_soma["Specificity"].mean(),
        metrics_combo["Specificity"].mean()
    ]
})

summary.to_csv(RESULTS_DIR/"model_comparison_summary.csv",index=False)

# =====================================================
# PANEL STABILITY
# =====================================================

panel_freq={}

for panel in outer_panels:
    for p in panel:
        panel_freq[p]=panel_freq.get(p,0)+1

panel_freq=(pd.Series(panel_freq)/OUTER_FOLDS)
panel_freq=panel_freq.sort_values(ascending=False)
panel_freq.to_csv(RESULTS_DIR/"outer_panel_frequency.csv")

print("\n=======================")
print(metrics_soma.mean(numeric_only=True))

metrics_soma.mean(numeric_only=True).to_csv(RESULTS_DIR/"metrics_mean.csv")
metrics_combo.mean(numeric_only=True).to_csv(RESULTS_DIR/"metrics_mean_combo.csv")

# =====================================================
# OOF EXPORT
# =====================================================

oof_df=pd.DataFrame({
"SampleId":sample_ids,
"y_true":y,
"oof_soma":oof_soma,
"oof_combo":oof_combo
})

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
"Model":["SomaScan","Combo"],
"Mean_CN":[oof_soma[y==0].mean(),oof_combo[y==0].mean()],
"Mean_AD":[oof_soma[y==1].mean(),oof_combo[y==1].mean()]
})

score_summary.to_csv(RESULTS_DIR/"score_summary.csv",index=False)

# =====================================================
# DELONG
# =====================================================

from scipy.stats import norm

assert oof_soma.notna().all()
assert oof_combo.notna().all()

def compute_midrank(x):
    order=np.argsort(x)
    z=x[order]
    T=np.zeros(len(x))
    i=0
    while i<len(x):
        j=i
        while j<len(x) and z[j]==z[i]:
            j+=1
        T[i:j]=0.5*(i+j-1)+1
        i=j
    out=np.empty(len(x))
    out[order]=T
    return out

def fast_delong(predictions_sorted_transposed,label_1_count):
    m=label_1_count
    n=predictions_sorted_transposed.shape[1]-m
    k=predictions_sorted_transposed.shape[0]
    tx=np.empty([k,m+n])

    for r in range(k):
        tx[r]=compute_midrank(predictions_sorted_transposed[r])
    aucs=(tx[:,:m].mean(axis=1)-tx[:,m:].mean(axis=1))/n
    v01=(tx[:,:m]-tx[:,:m].mean(axis=1)[:,None])/n
    v10=(tx[:,m:]-tx[:,m:].mean(axis=1)[:,None])/n
    cov=np.cov(v01)/m+np.cov(v10)/n
    return aucs,cov

def delong_pairwise(y_true,score_A,score_B):
    y_true=np.asarray(y_true)
    score_A=np.asarray(score_A)
    score_B=np.asarray(score_B)
    pos=np.where(y_true==1)[0]
    neg=np.where(y_true==0)[0]
    pred=np.vstack([
        np.concatenate([score_A[pos],score_A[neg]]),
        np.concatenate([score_B[pos],score_B[neg]])
    ])

    auc_delong,cov=fast_delong(pred,len(pos))
    auc1=roc_auc_score(y_true,score_A)
    auc2=roc_auc_score(y_true,score_B)

    var=cov[0,0]+cov[1,1]-2*cov[0,1]
    if var<1e-12:raise ValueError("Invalid DeLong variance")
    z=(auc1-auc2)/np.sqrt(var)
    p=2*(1-norm.cdf(abs(z)))
    return auc1,auc2,z,p

comp=[]

for name,a,b in [
("Soma_vs_Combo",oof_soma,oof_combo)
]:
    auc1,auc2,z,p=delong_pairwise(y.values,a,b)

    comp.append({
    "Comparison":name,
    "AUC_1":auc1,
    "AUC_2":auc2,
    "Delta_AUC":auc2-auc1,
    "Z":z,
    "P":p
    })

delong=pd.DataFrame(comp)
delong.to_csv(RESULTS_DIR/"delong_results.csv",index=False)

print("\n=======================")
print("\nDeLong")
print(delong)

# =====================================================
# FINAL MODEL + PERMUTATION IMPORTANCE
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

    best=pd.DataFrame(best_params_combo)
    kernel=best["kernel"].mode()[0]
    C=best["C"].median()
    gamma=best["gamma"].median()
    X_final=pd.concat([
    df[panel_final].reset_index(drop=True),
    X_apoe.reset_index(drop=True)],axis=1)
    scaler=StandardScaler()
    feature_names=X_final.columns
    X_final=scaler.fit_transform(X_final)
    svm=SVC(
        kernel=kernel,
        C=C,
        gamma=gamma,
        probability=False,
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
        "Proteins":feature_names,
        "Permutation":perm.importances_mean,
        "Permutation_sd":perm.importances_std
    })

    perm_summary=(perm_summary.sort_values("Permutation",ascending=False))
    perm_summary.to_csv(FINAL_PANEL_DIR/f"permutation_importance_{tag}.csv",index=False)
