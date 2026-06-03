# ============================================================
# Nested CV
# DEP -> Correlation -> Elastic Net -> SVM
# Nested CV
# ============================================================

import os
import numpy as np
import pandas as pd

from scipy.stats import spearmanr

from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score,
    roc_curve,
    accuracy_score,
    recall_score,
    precision_score,
    f1_score,
    balanced_accuracy_score,
    confusion_matrix
)

from sklearn.linear_model import LogisticRegressionCV
from sklearn.svm import SVC
from sklearn.inspection import permutation_importance
from skopt import BayesSearchCV
from skopt.space import Real, Categorical

# -----------------------------------------------------
# REPRODUCIBILITY
# -----------------------------------------------------

RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

# -----------------------------------------------------
# PATHS
# -----------------------------------------------------

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA_DIR = rf"{PROJECT_DIR}\Data"
RESULTS_DIR = rf"{PROJECT_DIR}\Results_proteins_Matched"
os.makedirs(RESULTS_DIR, exist_ok=True)
DATA_FILE = rf"{DATA_DIR}\Matched_Output.csv"

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
    "Sex",
    "Age",
    "ApoE",
    "Education",
    "p.tau217",
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
    "SampleGroup_bin",
    "distance",
    "weights",
    "subclass"    
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

outer_cv = StratifiedKFold(
    n_splits=OUTER_FOLDS,
    shuffle=True,
    random_state=RANDOM_STATE
)

metrics=[]
selection_summary=[]
outer_panels=[]
oof=np.zeros(len(df))
best_params=[]
roc_store=[]
STABILITY_THRESHOLD=0.6

# =====================================================
# LOOP INNER CV → STABILITY
# =====================================================

for fold,(train_idx,test_idx) in enumerate(outer_cv.split(X,y),start=1):

    print(f"\nFold {fold}")

    X_outer_train=X.iloc[train_idx].copy()
    X_outer_test=X.iloc[test_idx].copy()

    y_outer_train=y.iloc[train_idx]
    y_outer_test=y.iloc[test_idx]

    keep=X.columns.values

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
            solver="saga",
            Cs=np.logspace(-4,-1,30),
            l1_ratios=[0.90,0.95,1.0],
            cv=INNER_FOLDS,
            max_iter=MAX_ITER,
            tol=TOL,
            use_legacy_attributes=False,
            n_jobs=-1
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

    panel=freq[freq>=STABILITY_THRESHOLD].index.tolist()

    if len(panel)==0:
        panel=freq.sort_values(ascending=False).head(10).index.tolist()

    outer_panels.append(panel)  
    pd.Series(panel).to_csv(rf"{RESULTS_DIR}\panel_fold_{fold}.csv",index=False)
    n_selected=len(panel)   
    selection_summary.append({"Fold":fold,"N_features":n_selected})

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
        n_iter=30,
        cv=INNER_FOLDS,
        scoring="roc_auc",
        n_jobs=-1,
        random_state=RANDOM_STATE
    )

    search.fit(X_train,y_outer_train)
    best_params.append(search.best_params_)
    y_prob=search.best_estimator_.decision_function(X_test)
    oof[test_idx]=(y_prob-y_prob.mean())/(y_prob.std()+1e-8)
    y_pred=(y_prob>=0).astype(int)
    tn,fp,fn,tp=confusion_matrix(y_outer_test,y_pred).ravel()
    
    auc=roc_auc_score(y_outer_test,y_prob)
    fpr,tpr,thr=roc_curve(y_outer_test,y_prob)
    roc_store.append(pd.DataFrame({
    "Fold":fold,"Model":"Proteins","FPR":fpr,"TPR":tpr,"Threshold":thr,"AUC":auc}))
       
    metrics.append({
    "Fold":fold,
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
metrics_df.to_csv(rf"{RESULTS_DIR}\metrics_nestedCV.csv",index=False)
pd.DataFrame(selection_summary).to_csv(rf"{RESULTS_DIR}\feature_counts.csv",index=False)

print("\nOOF AUC")
print(roc_auc_score(y,oof))
print("\nSaved.")

# ==========================================
# PANEL STABILITY
# ==========================================

panel_freq={}

for panel in outer_panels:
    for p in panel:
        panel_freq[p]=panel_freq.get(p,0)+1

panel_freq=(pd.Series(panel_freq)/OUTER_FOLDS)
panel_freq=panel_freq.sort_values(ascending=False)
panel_freq.to_csv(rf"{RESULTS_DIR}\outer_panel_frequency.csv")

print("\n=======================")
print(metrics_df.mean(numeric_only=True))

# =====================================================
# OOF EXPORT
# =====================================================

oof_df=pd.DataFrame({
"SampleId":df.index,
"y_true":y,
"oof":oof
})

oof_df.to_csv(rf"{RESULTS_DIR}\oof_predictions.csv",index=False)

# =====================================================
# ROC
# =====================================================

roc_df=pd.concat(roc_store,ignore_index=True)
roc_df.to_csv(rf"{RESULTS_DIR}\roc_curves.csv",index=False)

# =====================================================
# SCORE SUMMARY
# =====================================================

score_summary=pd.DataFrame({
"Model":["SomaScan"],"Mean_CN":[oof[y==0].mean()],"Mean_AD":[oof[y==1].mean()]})

score_summary.to_csv(rf"{RESULTS_DIR}\score_summary.csv",index=False)

# =====================================================
# FINAL MODEL + PERMUTATION
# =====================================================

panel_final=(panel_freq[panel_freq>=STABILITY_THRESHOLD].index.tolist())

if len(panel_final)==0:
    panel_final=(panel_freq.head(10).index.tolist())

print("\n=======================")
print("\nFinal panel")
print(panel_final)

pd.Series(panel_final).to_csv(rf"{RESULTS_DIR}\panel_final.csv",index=False)

best=pd.DataFrame(best_params)
best.to_csv(rf"{RESULTS_DIR}\best_params_outer.csv",index=False) #NUEVO
kernel=best["kernel"].mode()[0]
C=best["C"].median()
gamma=best["gamma"].median()

print("\nFinal SVM")
print(kernel)
print(C)
print(gamma)

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

perm_summary=(perm_summary.sort_values("Permutation",ascending=False))
perm_summary.to_csv(rf"{RESULTS_DIR}\permutation_importance.csv",index=False)

print("\nPermutation")
print(perm_summary)

