###############################################################################
# ReDLat plasma proteomics — strict machine-learning workflow
# 13. Evaluate the fixed panel with p-tau217
# Requires: private master matrix and prespecified seven-protein panel
# Produces: strict OOF logistic comparisons and DeLong tests
# Data policy: participant-level files remain under the private result root.
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


import numpy as np
import pandas as pd
from scipy.stats import chi2
import statsmodels.api as sm
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from redlat_ml.config import load_config, require_files
from redlat_ml.strict_cv import StrictCVSettings, delong_pairwise, resolve_column, safe_n_splits

CONFIG = load_config(__file__); SETTINGS = StrictCVSettings.from_env()
RESULTS_DIR = CONFIG.private_root / "fixed_panel_ptau"; RESULTS_DIR.mkdir(parents=True, exist_ok=True)
require_files([(CONFIG.master_file, "ML master matrix")])
PANEL = ["SPC25", "CPLX2", "TCP11L1", "ACHE", "ODC1", "SPON1", "RTN4RL1"]
df = pd.read_csv(CONFIG.master_file, low_memory=False); df = df[df["SampleGroup"].isin(["CN", "AD"])].copy()
ptau = resolve_column(df, ("p.tau217", "p-tau217", "pTau217", "ptau217")); df[ptau] = pd.to_numeric(df[ptau], errors="coerce")
for protein in PANEL:
    if protein not in df.columns: raise KeyError(f"Fixed-panel protein not found: {protein}")
flow = pd.DataFrame({"Stage": ["CN/AD cohort", "p-tau217 available"], "N": [len(df), int(df[ptau].notna().sum())], "CN": [int((df.SampleGroup == 'CN').sum()), int(((df.SampleGroup == 'CN') & df[ptau].notna()).sum())], "AD": [int((df.SampleGroup == 'AD').sum()), int(((df.SampleGroup == 'AD') & df[ptau].notna()).sum())], "Protein_complete_case_required": [False, False]})
df = df[df[ptau].notna()].reset_index(drop=True); y = df.SampleGroup.map({"CN": 0, "AD": 1}).astype(int)
models = {"Panel": PANEL, "pTau217": [ptau], "Panel+pTau217": PANEL + [ptau]}
outer = StratifiedKFold(n_splits=safe_n_splits(y, 5), shuffle=True, random_state=SETTINGS.random_state)
oof = pd.DataFrame(index=df.index); oof["y_true"] = y; roc_rows = []; auc_rows = []
for model, features in models.items():
    probabilities = np.full(len(df), np.nan); fold_aucs = []
    for fold, (train, test) in enumerate(outer.split(df, y), start=1):
        inner = StratifiedKFold(n_splits=safe_n_splits(y.iloc[train], SETTINGS.inner_folds), shuffle=True, random_state=SETTINGS.random_state + fold)
        pipe = Pipeline([("imputer", SimpleImputer(strategy="median")), ("scaler", StandardScaler()), ("classifier", LogisticRegression(max_iter=SETTINGS.max_iter, random_state=SETTINGS.random_state))])
        search = GridSearchCV(pipe, {"classifier__C": np.logspace(-3, 3, 13)}, scoring="roc_auc", cv=inner, n_jobs=SETTINGS.n_jobs, refit=True)
        search.fit(df.iloc[train][features], y.iloc[train]); probability = search.predict_proba(df.iloc[test][features])[:, 1]; probabilities[test] = probability
        fpr, tpr, threshold = roc_curve(y.iloc[test], probability); fold_auc = roc_auc_score(y.iloc[test], probability); fold_aucs.append(fold_auc)
        roc_rows.append(pd.DataFrame({"Model": model, "Fold": fold, "FPR": fpr, "TPR": tpr, "Threshold": threshold, "AUC": fold_auc}))
    oof[model] = probabilities; auc_rows.append({"Model": model, "AUC": roc_auc_score(y, probabilities), "Mean_fold_AUC": np.mean(fold_aucs), "SD_fold_AUC": np.std(fold_aucs, ddof=1), "N": len(df)})
pd.DataFrame(auc_rows).to_csv(RESULTS_DIR / "auc_summary.csv", index=False); pd.concat(roc_rows, ignore_index=True).to_csv(RESULTS_DIR / "roc_fold_curves.csv", index=False); oof.to_csv(RESULTS_DIR / "oof_predictions.csv", index=False); flow.to_csv(RESULTS_DIR / "cohort_summary.csv", index=False)
comparisons = []
for label, a, b in (("Panel_vs_pTau", "Panel", "pTau217"), ("Panel_vs_Combo", "Panel", "Panel+pTau217"), ("pTau_vs_Combo", "pTau217", "Panel+pTau217")):
    comparisons.append({"Comparison": label, **delong_pairwise(y, oof[a], oof[b])})
pd.DataFrame(comparisons).to_csv(RESULTS_DIR / "delong_results.csv", index=False)

# Descriptive complete-case likelihood models; these do not contribute to OOF AUC.
stats_rows, fits = [], {}
for model, features in models.items():
    work = df[features].apply(pd.to_numeric, errors="coerce").dropna(); y_complete = y.loc[work.index]
    scaled = StandardScaler().fit_transform(work); design = sm.add_constant(pd.DataFrame(scaled, columns=features, index=work.index), has_constant="add")
    fit = sm.Logit(y_complete, design).fit(disp=False); fits[model] = fit
    stats_rows.append({"Model": model, "McFadden_R2": fit.prsquared, "AIC": fit.aic, "BIC": fit.bic, "LogLik": fit.llf, "N": len(work), "Role": "descriptive_full_cohort_model"})
pd.DataFrame(stats_rows).to_csv(RESULTS_DIR / "model_statistics.csv", index=False)
lr1 = 2 * (fits["Panel+pTau217"].llf - fits["Panel"].llf); df1 = fits["Panel+pTau217"].df_model - fits["Panel"].df_model
lr2 = 2 * (fits["Panel+pTau217"].llf - fits["pTau217"].llf); df2 = fits["Panel+pTau217"].df_model - fits["pTau217"].df_model
pd.DataFrame({"Comparison": ["Panel_vs_Combo", "pTau_vs_Combo"], "LR": [lr1, lr2], "df": [df1, df2], "P": [chi2.sf(lr1, df1), chi2.sf(lr2, df2)], "Role": ["descriptive"] * 2}).to_csv(RESULTS_DIR / "likelihood_ratio_tests.csv", index=False)
df[PANEL + [ptau]].corr(method="spearman").to_csv(RESULTS_DIR / "correlation_matrix.csv")
print(pd.DataFrame(auc_rows))
