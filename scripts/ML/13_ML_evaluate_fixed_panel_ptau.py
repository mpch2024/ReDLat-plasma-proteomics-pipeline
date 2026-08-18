###############################################################################
# ReDLat plasma proteomics — machine-learning workflow
# 13. Evaluate the fixed panel with p-tau217
# Requires: private master matrix and the primary seven-protein panel
# Produces: cross-validated logistic comparisons and DeLong tests
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

import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score, roc_curve
import statsmodels.api as sm
from scipy.stats import chi2

# =====================================================
# SETTINGS
# =====================================================

RANDOM_STATE = 42
N_SPLITS = 5

PROJECT_DIR = CONFIG.project_root
DATA_FILE = CONFIG.master_file
RESULTS_DIR = CONFIG.private_root / "fixed_panel_ptau"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
TARGET = "SampleGroup"
PTAU = "p.tau217"
FINAL_PANEL = ["SPC25","CPLX2","TCP11L1","ACHE","ODC1","SPON1","RTN4RL1"]

# =====================================================
# LOAD DATA
# =====================================================

df = pd.read_csv(DATA_FILE)
df = df[df[TARGET].isin(["CN", "AD"])].copy()
COMMON_COLS = FINAL_PANEL + [PTAU]
N0 = len(df)
df = df.dropna(subset=COMMON_COLS).copy()
N1 = len(df)

print("\n========================")
print(f"Original N = {N0}")
print(f"Final N    = {N1}")
print(f"Removed    = {N0-N1}")
print("========================\n")

pd.DataFrame({
    "Original_N":[N0],
    "Final_N":[N1],
    "Removed":[N0-N1]
}).to_csv(
    RESULTS_DIR/"cohort_summary.csv",
    index=False
)

y = df[TARGET].map({"CN":0,"AD":1}).values

# =====================================================
# MODELS
# =====================================================

MODELS = {
    "Panel": FINAL_PANEL,
    "pTau217": [PTAU],
    "Panel+pTau217": FINAL_PANEL + [PTAU]
}

# =====================================================
# 5-FOLD CV
# =====================================================
roc_fold_rows = []

cv = StratifiedKFold(
    n_splits=N_SPLITS,
    shuffle=True,
    random_state=RANDOM_STATE
)

auc_rows = []
roc_rows = []

oof_df = pd.DataFrame(index=df.index)
oof_df["y_true"] = y

for model_name, features in MODELS.items():

    print(f"\nRunning: {model_name}")

    oof_prob = np.zeros(len(df))

    for fold, (train_idx, test_idx) in enumerate(
        cv.split(df, y),
        start=1
    ):

        X_train = df.iloc[train_idx][features]
        X_test = df.iloc[test_idx][features]

        pipe = Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(
                penalty="l2",
                C=1.0,
                max_iter=10000,
                random_state=RANDOM_STATE
            ))
        ])

        pipe.fit(X_train,y[train_idx])
        prob = pipe.predict_proba(X_test)[:,1]
        fpr_fold, tpr_fold, _ = roc_curve(y[test_idx],prob)

        roc_fold_rows.append(
            pd.DataFrame({
                "Model": model_name,
                "Fold": fold,
                "FPR": fpr_fold,
                "TPR": tpr_fold
            })
        )

        oof_prob[test_idx] = prob
    auc_val = roc_auc_score(y,oof_prob)

    print(f"AUC = {auc_val:.3f}")

    fpr, tpr, _ = roc_curve(y,oof_prob)
    auc_rows.append({"Model":model_name,"AUC":auc_val})
    roc_tmp = pd.DataFrame({"Model":model_name,"FPR":fpr,"TPR":tpr})
    roc_rows.append(roc_tmp)
    oof_df[model_name] = oof_prob

# =====================================================
# EXPORT CV RESULTS
# =====================================================

pd.DataFrame(auc_rows).to_csv(RESULTS_DIR/"auc_summary.csv",index=False)
pd.concat(roc_fold_rows,ignore_index=True).to_csv(RESULTS_DIR/"roc_fold_curves.csv",index=False)
oof_df.to_csv(RESULTS_DIR/"oof_predictions.csv",index=False)

# =====================================================
# FULL LOGISTIC MODELS
# =====================================================

stats_rows = []
fitted_models = {}
for model_name, features in MODELS.items():
    X = df[features].copy()
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    X_scaled = pd.DataFrame(X_scaled,columns=features,index=df.index)
    X_scaled = sm.add_constant(X_scaled,has_constant="add")
    fit = sm.Logit(y,X_scaled).fit(disp=False)
    fitted_models[model_name] = fit
    stats_rows.append({
        "Model": model_name,
        "McFadden_R2": fit.prsquared,
        "AIC": fit.aic,
        "BIC": fit.bic,
        "LogLik": fit.llf,
        "N": len(df)
    })

    ci = fit.conf_int()

    coef_table = pd.DataFrame({
        "Variable": fit.params.index,
        "Beta": fit.params.values,
        "OR": np.exp(fit.params.values),
        "CI_low": np.exp(ci.iloc[:,0]),
        "CI_high": np.exp(ci.iloc[:,1]),
        "P": fit.pvalues.values
    })

    safe_name = (model_name.replace("+","_plus_").replace(" ","_"))
    coef_table.to_csv(RESULTS_DIR /f"coefficients_{safe_name}.csv",index=False)

# =====================================================
# MODEL STATISTICS
# =====================================================

pd.DataFrame(stats_rows).to_csv(RESULTS_DIR/"model_statistics.csv",index=False)

# =====================================================
# LIKELIHOOD RATIO TESTS
# =====================================================

panel_fit = fitted_models["Panel"]
ptau_fit = fitted_models["pTau217"]
combo_fit = fitted_models["Panel+pTau217"]

# ---------------------------
# pTau adds to panel
# ---------------------------

lr1 = 2 * (combo_fit.llf - panel_fit.llf)
df1 = (combo_fit.df_model - panel_fit.df_model)
p1 = chi2.sf(lr1,df1)

# ---------------------------
# panel adds to pTau
# ---------------------------

lr2 = 2 * (combo_fit.llf - ptau_fit.llf)
df2 = (combo_fit.df_model - ptau_fit.df_model)
p2 = chi2.sf(lr2,df2)

pd.DataFrame({
    "Comparison":["Panel_vs_Combo","pTau_vs_Combo"],"LR":[lr1,lr2],"df":[df1,df2],"P":[p1,p2]
}).to_csv(RESULTS_DIR/"likelihood_ratio_tests.csv",index=False)

# =====================================================
# CORRELATIONS
# =====================================================

corr_df = df[FINAL_PANEL + [PTAU]].corr(method="spearman")
corr_df.to_csv(RESULTS_DIR/"correlation_matrix.csv")

print("\nAnalysis completed.")

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

oof_panel = oof_df["Panel"]
oof_ptau = oof_df["pTau217"]
oof_combo = oof_df["Panel+pTau217"]

y_true = oof_df["y_true"]

comp=[]

for name,a,b in [
("Panel_vs_pTau",oof_panel,oof_ptau),
("Panel_vs_Combo",oof_panel,oof_combo),
("pTau_vs_Combo",oof_ptau,oof_combo)
]:
    auc1,auc2,z,p=delong_pairwise(y_true,a,b)

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