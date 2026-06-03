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
RESULTS_DIR = rf"{PROJECT_DIR}\Results_proteins_with_ApoE"
os.makedirs(RESULTS_DIR, exist_ok=True)
DATA_FILE = rf"{DATA_DIR}\somascan_filter_DEP_final_FDR_metadata.csv"

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

df=pd.read_csv(DATA_FILE)

META_COLS=[
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
    "T.ADLQ"
]

df[ID_COL]=df[ID_COL].astype(str)

df=df.set_index(ID_COL)

df=df[df[TARGET].isin(["CN","AD"])].copy()

protein_cols=[
    c
    for c in df.columns
    if c not in META_COLS
]

# =====================================================
# COMMON COHORT
# proteins + ApoE only
# =====================================================

N0=len(df)
COMMON_COLS=protein_cols+["ApoE"]
df=df.dropna(subset=COMMON_COLS).reset_index()
N1=len(df)

pd.DataFrame({
    "Original_N":[N0],
    "Final_N":[N1],
    "Removed":[N0-N1]
}).to_csv(rf"{RESULTS_DIR}\selected_subjects.csv",index=False)

# =====================================================
# APOE ε4 carrier
# =====================================================

df["ApoE"]=(df["ApoE"].astype(str).str.contains("4",na=False).astype(int))

df[ID_COL]=df[ID_COL].astype(str)

sample_ids=df[ID_COL].values

y=df[TARGET].map({"CN":0,"AD":1})

# =====================================================
# MATRICES
# =====================================================

X_soma=df[protein_cols]
X_apoe=df[["ApoE"]]

print("\nSoma:",X_soma.shape)
print("\nApoE:",X_apoe.shape)

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

outer_cv=StratifiedKFold(
    n_splits=OUTER_FOLDS,
    shuffle=True,
    random_state=RANDOM_STATE
)

metrics_soma=[]
metrics_combo=[]
best_params_soma=[]
best_params_combo=[]
outer_panels=[]
roc_store=[]

oof_soma=np.zeros(len(df))
oof_combo=np.zeros(len(df))

STABILITY_THRESHOLD=0.8

print("\nApoE overall")

print(df.groupby(TARGET)["ApoE"].mean())

# =====================================================
# LOOP
# =====================================================

for fold,(train_idx,test_idx) in enumerate(
    outer_cv.split(X_soma,y),start=1):

    print(f"\nFold {fold}")

    # =================================================
    # SOMASCAN
    # =================================================

    X_outer_train=X_soma.iloc[train_idx].copy()
    X_outer_test=X_soma.iloc[test_idx].copy()
    
    y_outer_train=y.iloc[train_idx]
    y_outer_test=y.iloc[test_idx]

    keep=X_soma.columns.values

    inner_selected=[]

    inner_cv=StratifiedKFold(
        n_splits=INNER_FOLDS,
        shuffle=True,
        random_state=RANDOM_STATE
    )

    for itr,ival in inner_cv.split(
        X_outer_train,
        y_outer_train
    ):

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

    # CONSENSUS PANEL
    freq={}

    for s in inner_selected:
        for p in s:
            freq[p]=(freq.get(p,0)+1)

    freq=(pd.Series(freq)/INNER_FOLDS)

    panel=freq[freq>=STABILITY_THRESHOLD].index.tolist()

    if len(panel)==0:
        panel=(freq.sort_values(ascending=False).head(10).index.tolist())
    
    outer_panels.append(panel)

    pd.Series(panel).to_csv(rf"{RESULTS_DIR}\panel_fold_{fold}.csv",index=False)

    # SOMASCAN
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
    best_params_soma.append(search.best_params_)
    y_prob=search.best_estimator_.decision_function(X_test)           
    oof_soma[test_idx]=(y_prob-y_prob.mean())/(y_prob.std()+1e-8)
    
    y_pred=(y_prob>=0).astype(int)
    tn,fp,fn,tp=confusion_matrix(y_outer_test,y_pred).ravel()
    
    auc=roc_auc_score(y_outer_test,y_prob)
    fpr,tpr,thr=roc_curve(y_outer_test,y_prob)
    roc_store.append(pd.DataFrame({
    "Fold":fold,"Model":"Proteins","FPR":fpr,"TPR":tpr,"Threshold":thr,"AUC":auc}))   
    
    metrics_soma.append({
    "Fold":fold,
    "AUC":roc_auc_score(y_outer_test,y_prob),
    "Accuracy":accuracy_score(y_outer_test,y_pred),
    "Sensitivity":recall_score(y_outer_test,y_pred),
    "Specificity":tn/(tn+fp),
    "Precision":precision_score(y_outer_test,y_pred,zero_division=np.nan),
    "F1":f1_score(y_outer_test,y_pred,zero_division=np.nan)
    })

    print(f"Soma AUC={metrics_soma[-1]['AUC']:.3f}")
    
    # =================================================
    # SOMASCAN + ApoE
    # =================================================
    
    X_train_combo=X_outer_train[panel].copy()
    X_train_combo["ApoE"]=(X_apoe.iloc[train_idx].values)
    
    X_test_combo=X_outer_test[panel].copy()
    X_test_combo["ApoE"]=(X_apoe.iloc[test_idx].values)
    
    scaler=StandardScaler()      
    
    X_train_combo=scaler.fit_transform(X_train_combo)
    X_test_combo=scaler.transform(X_test_combo)

    search=BayesSearchCV(
    estimator=SVC(probability=False,random_state=RANDOM_STATE),
    search_spaces=svm_param_space,
    n_iter=30,
    cv=INNER_FOLDS,
    scoring="roc_auc",
    n_jobs=-1,
    random_state=RANDOM_STATE
    )

    search.fit(X_train_combo,y_outer_train)
    best_params_combo.append(search.best_params_)    
    y_prob=search.best_estimator_.decision_function(X_test_combo)        
    oof_combo[test_idx]=(y_prob-y_prob.mean())/(y_prob.std()+1e-8)      
    y_pred=(y_prob>=0).astype(int)
    tn,fp,fn,tp=confusion_matrix(y_outer_test,y_pred).ravel()
    
    auc=roc_auc_score(y_outer_test,y_prob)
    fpr,tpr,thr=roc_curve(y_outer_test,y_prob)
    roc_store.append(pd.DataFrame({
    "Fold":fold,"Model":"Proteins + ApoE","FPR":fpr,"TPR":tpr,"Threshold":thr,"AUC":auc}))

    metrics_combo.append({
    "Fold":fold,
    "AUC":auc,
    "Accuracy":accuracy_score(y_outer_test,y_pred),
    "Sensitivity":recall_score(y_outer_test,y_pred),
    "Specificity":tn/(tn+fp),
    "Precision":precision_score(y_outer_test,y_pred,zero_division=np.nan),
    "F1":f1_score(y_outer_test,y_pred,zero_division=np.nan)
    })

    print(f"Combo AUC={metrics_combo[-1]['AUC']:.3f}")
    
# =====================================================
# SAVE METRICS
# =====================================================

metrics_soma=pd.DataFrame(metrics_soma)
metrics_combo=pd.DataFrame(metrics_combo)

metrics_soma.to_csv(rf"{RESULTS_DIR}\metrics_soma.csv",index=False)
metrics_combo.to_csv(rf"{RESULTS_DIR}\metrics_combo.csv",index=False)

summary=pd.concat([
metrics_soma.mean(numeric_only=True).to_frame().T.assign(Model="SomaScan"),
metrics_combo.mean(numeric_only=True).to_frame().T.assign(Model="Combo")
])

summary=summary[[
"Model",
"AUC",
"Accuracy",
"Sensitivity",
"Specificity",
"Precision",
"F1"
]]

summary.columns=[
"Model",
"AUC_mean",
"Accuracy_mean",
"Sensitivity_mean",
"Specificity_mean",
"Precision_mean",
"F1_mean"
]

summary.to_csv(rf"{RESULTS_DIR}\metrics_summary.csv",index=False)

print("\n=======================")
print(summary)

print("\nOOF SOMASCAN")
print(roc_auc_score(y,oof_soma))
print(pd.Series(oof_soma).describe())
print("\nOOF COMBO")
print(roc_auc_score(y,oof_combo))
print(pd.Series(oof_combo).describe())


# =====================================================
# PANEL STABILITY
# =====================================================

panel_freq={}

for panel in outer_panels:
    for p in panel:
        panel_freq[p]=panel_freq.get(p,0)+1

panel_freq=(pd.Series(panel_freq)/OUTER_FOLDS)
panel_freq=panel_freq.sort_values(ascending=False)

panel_freq.to_csv(rf"{RESULTS_DIR}\outer_panel_frequency.csv")

# =====================================================
# OOF EXPORT
# =====================================================

oof_df=pd.DataFrame({
"SampleId":sample_ids,
"y_true":y,
"oof_soma":oof_soma,
"oof_combo":oof_combo
})

oof_df.to_csv(rf"{RESULTS_DIR}\oof_predictions.csv",index=False)

# =====================================================
# ROC
# =====================================================

roc_df=pd.concat(roc_store,ignore_index=True)
roc_df.to_csv(rf"{RESULTS_DIR}\roc_curves.csv",index=False)

#fpr_soma,tpr_soma,thr_soma=roc_curve(y,oof_soma)
#fpr_combo,tpr_combo,thr_combo=roc_curve(y,oof_combo)
#auc_soma=roc_auc_score(y,oof_soma)
#auc_combo=roc_auc_score(y,oof_combo)
#roc_df=pd.concat([
#pd.DataFrame({"Model":"SomaScan","FPR":fpr_soma,"TPR":tpr_soma,"Threshold":thr_soma,"AUC":auc_soma}),
#pd.DataFrame({"Model":"Combo","FPR":fpr_combo,"TPR":tpr_combo,"Threshold":thr_combo,"AUC":auc_combo})])
#roc_df.to_csv(rf"{RESULTS_DIR}\roc_curves.csv",index=False)

# =====================================================
# SCORE SUMMARY
# =====================================================

score_summary=pd.DataFrame({
"Model":["SomaScan","Combo"],
"Mean_CN":[oof_soma[y==0].mean(),oof_combo[y==0].mean()],
"Mean_AD":[oof_soma[y==1].mean(),oof_combo[y==1].mean()]
})

score_summary.to_csv(rf"{RESULTS_DIR}\score_summary.csv",index=False)

# =====================================================
# DELONG
# =====================================================

from scipy.stats import norm

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
    y_true=np.array(y_true)
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

auc1,auc2,z,p=delong_pairwise(y.values,oof_soma,oof_combo)

comp.append({
"Comparison":"Soma_vs_ApoECombo",
"AUC_1":auc1,
"AUC_2":auc2,
"Delta_AUC":auc2-auc1,
"Z":z,
"P":p
})

delong=pd.DataFrame(comp)
delong.to_csv(rf"{RESULTS_DIR}\delong_results.csv",index=False)

print("\n=======================")
print("\nDeLong")
print(delong)

# =====================================================
# FINAL MODEL + PERMUTATION IMPORTANCE
# =====================================================

panel_final=(panel_freq[panel_freq>=STABILITY_THRESHOLD].index.tolist())

if len(panel_final)==0:
    panel_final=(panel_freq.head(10).index.tolist())

pd.Series(panel_final).to_csv(rf"{RESULTS_DIR}\panel_final.csv",index=False)
best=pd.DataFrame(best_params_combo)
best.to_csv(rf"{RESULTS_DIR}\best_params_outer_combo.csv",index=False)

kernel=best["kernel"].mode()[0]
C=best["C"].median()
gamma=best["gamma"].median()

X_final=df[panel_final].copy()
X_final["ApoE"]=(df["ApoE"].values)
feature_names=(X_final.columns)

scaler=StandardScaler()

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
"Protein":feature_names,
"Permutation_mean":perm.importances_mean,
"Permutation_sd":perm.importances_std
})

perm_summary=(perm_summary.sort_values("Permutation_mean",ascending=False))
perm_summary.to_csv(rf"{RESULTS_DIR}\permutation_importance.csv",index=False)

print("\nTop Permutation")
print(perm_summary)