###############################################################################
# ReDLat plasma proteomics — machine-learning workflow
# 14. Model clinical alignment
# Requires: private master matrix and the fixed seven-protein panel
# Produces: regression metrics and coefficient tables
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

import pandas as pd
import numpy as np
import os
from pathlib import Path
from sklearn.model_selection import KFold
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import scipy.stats
import math
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import warnings
warnings.filterwarnings("ignore")


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
RESULTS_DIR = CONFIG.private_root / "clinical_regression"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
DATA_FILE = CONFIG.master_file

# -----------------------------------------------------
# LOAD
# -----------------------------------------------------

df_raw = pd.read_csv(DATA_FILE)
df_raw = df_raw[df_raw['SampleGroup'].isin(['CN', 'AD'])].copy()

# -----------------------------------------------------
# Variables
# -----------------------------------------------------

protein_cols =["SPC25","CPLX2","TCP11L1","ACHE","ODC1","SPON1","RTN4RL1"]
predictors = protein_cols
outcomes = {
    'mmse_total': 'Global Cognition',
    'cog_craft_verb_delayed': 'Memory',
    'udsfaq_total': 'Functionality',
    'T.ADLQ': 'Deterioration',
    'cdr_boxscore': 'Dementia Severity'
}

# -----------------------------------------------------
# Prepare predictors
# -----------------------------------------------------

X = df_raw[predictors].copy()

print(X.columns)

# Reverse scales
df_raw['udsfaq_total'] = df_raw['udsfaq_total'].max() - df_raw['udsfaq_total']
df_raw['T.ADLQ'] = df_raw['T.ADLQ'].max() - df_raw['T.ADLQ']
df_raw['cdr_boxscore'] = df_raw['cdr_boxscore'].max() - df_raw['cdr_boxscore']

# -----------------------------------------------------
# OLS utilities
# -----------------------------------------------------
def ols_t_pvalues(X_df, y_series):
    Xc = sm.add_constant(X_df, has_constant='add')
    model = sm.OLS(y_series, Xc).fit()
    tvals = pd.Series(index=['_intercept'] + list(X_df.columns), dtype=float)
    pvals = pd.Series(index=['_intercept'] + list(X_df.columns), dtype=float)
    tvals['_intercept'] = model.tvalues['const']
    pvals['_intercept'] = model.pvalues['const']
    for col in X_df.columns:
        tvals[col] = model.tvalues[col]
        pvals[col] = model.pvalues[col]

    return tvals, pvals

# -----------------------------------------------------
# Regression with outer CV + OLS
# -----------------------------------------------------
def Regression_Linear_outer_with_OLS(X, y, n_outer=10, random_state=42):

    p = X.shape[1]
    coef_array = np.zeros((p + 1, n_outer))
    r2_list, mse_list, rmse_list = [], [], []
    y_pred_all, y_test_all = [], []
    outer_cv = KFold(n_splits=n_outer, shuffle=True, random_state=random_state)

    for k, (tr, te) in enumerate(outer_cv.split(X, y)):
        X_tr, X_te = X.iloc[tr], X.iloc[te]
        y_tr, y_te = y.iloc[tr], y.iloc[te]

        scaler = MinMaxScaler(feature_range=(0.05, 0.95))
        X_tr_s = scaler.fit_transform(X_tr)
        X_te_s = scaler.transform(X_te)

        model = LinearRegression()
        model.fit(X_tr_s, y_tr)
        y_hat = model.predict(X_te_s)

        coef_array[0, k] = model.intercept_
        coef_array[1:, k] = model.coef_

        r2_list.append(r2_score(y_te, y_hat))
        mse_list.append(mean_squared_error(y_te, y_hat))
        rmse_list.append(math.sqrt(mean_squared_error(y_te, y_hat)))

        y_pred_all.extend(y_hat)
        y_test_all.extend(y_te.values)

    y_pred_all = np.array(y_pred_all)
    y_test_all = np.array(y_test_all)

    n = len(y_test_all)
    r2 = r2_score(y_test_all, y_pred_all)
    r2_adj = 1 - (1 - r2) * (n - 1) / (n - p - 1)

    mse = mean_squared_error(y_test_all, y_pred_all)
    rmse = math.sqrt(mse)

    F = (r2 / p) / ((1 - r2) / (n - p - 1))
    F_p = scipy.stats.f.sf(F, p, n - p - 1)
    F2 = r2 / (1 - r2)

    coef_mean = coef_array.mean(axis=1)
    coef_std = coef_array.std(axis=1)

    # ---- OLS FINAL
    scaler_full = StandardScaler()
    X_scaled = pd.DataFrame(scaler_full.fit_transform(X), columns=X.columns)
    y_aligned = y.reset_index(drop=True)

    tvals, pvals = ols_t_pvalues(X_scaled, y_aligned)

    coef_df = pd.DataFrame(
        index=['_intercept'] + list(X.columns),
        data={
            'Estimate mean': coef_mean,
            'Estimate std': coef_std,
            't value': tvals,
            'p value': pvals
        }
    )

    coef_df.loc['_intercept', 'R2'] = r2
    coef_df.loc['_intercept', 'R2 adj'] = r2_adj
    coef_df.loc['_intercept', 'F2'] = F2
    coef_df.loc['_intercept', 'F'] = F
    coef_df.loc['_intercept', 'F-p_value'] = F_p
    coef_df.loc['_intercept', 'rmse'] = rmse
    coef_df.loc['_intercept', 'outcome var'] = np.var(y, ddof=1)

    return coef_df


# -----------------------------------------------------
# Run regressions
# -----------------------------------------------------
coef_resultados = {}

for out, label in outcomes.items():
    print(f"\n== {label} ==")

    df = df_raw[predictors + [out]].dropna()

    # Separar predictores
    X = df[predictors].copy()

    # Convertir todo a numérico
    X = X.astype(float)

     # Outcome
    y = df[out].astype(float)

    # Correr regresión
    coef_df = Regression_Linear_outer_with_OLS(X, y)
    pd.DataFrame(coef_df).to_csv(RESULTS_DIR/f"Regression_CN_AD_{label.replace(' ', '')}.csv",index=False)
    
    # Store for plotting
    coef_resultados[label] = coef_df

# -----------------------------------------------------
# 9. Plot betas values
# -----------------------------------------------------
import math

proteins_up = ["SPC25","CPLX2","TCP11L1","ACHE","ODC1","SPON1"]
proteins_down = ["RTN4RL1"]
proteins = proteins_up + proteins_down

# Plot size
fig_width_cm = 18
fig_height_cm = 5.4
fig_size_inch = (fig_width_cm / 2.54, fig_height_cm / 2.54)

fig, axes = plt.subplots(
    1, len(coef_resultados),
    figsize=fig_size_inch,
    dpi=300)

if len(coef_resultados) == 1:
    axes = [axes]

for ax, (plot_title, coef_df) in zip(axes, coef_resultados.items()):

    # ONLY protein betas
    df_betas = coef_df.loc[proteins,['Estimate mean', 'p value']].copy()
    df_betas = df_betas.sort_values(by='Estimate mean')
    df_betas['label'] = df_betas.index

    # Colors
    colors = []
    for protein, pval in zip(df_betas.index, df_betas['p value']):
        if protein in proteins_up:
            colors.append('#d90f0f' if pval < 0.05 else '#edb9b9')
        else:
            colors.append('#1f77b4' if pval < 0.05 else '#aec7e8')

    sns.barplot(
        data=df_betas,
        y='label',
        x='Estimate mean',
        palette=colors,
        ax=ax
    )

    max_beta = df_betas['Estimate mean'].abs().max()
    ax.set_xlim(-max_beta * 1.1, max_beta * 1.1)
    ax.axvline(0, color='black', linestyle='--', linewidth=0.5)

    ax.set_title(plot_title, fontsize=7, pad=18)
    ax.set_xlabel('Beta value', fontsize=6)
    ax.set_ylabel('')
    ax.tick_params(axis='both', which='major',
               labelsize=5, pad=1,length=2, width=0.5)

    # Metrics
    r2 = coef_df.loc['_intercept', 'R2']
    f2 = coef_df.loc['_intercept', 'F2']
    fval = coef_df.loc['_intercept', 'F']
    fp = coef_df.loc['_intercept', 'F-p_value']

    metrics_line1 = f"R² = {r2:.2f}, F² = {f2:.2f}"
    metrics_line2 = (f"F = {fval:.2f}, "f"{'p < 1.0e-15' if fp < 1.0e-15 else f'p = {fp:.1e}'}")
    ax.text(0.5, 1.12, metrics_line1, fontsize=6, ha='center', va='center', transform=ax.transAxes)
    ax.text(0.5, 1.06, metrics_line2, fontsize=6, ha='center', va='center', transform=ax.transAxes)

plt.tight_layout()
plt.savefig(RESULTS_DIR/"Betas_Cognitive_Domains.pdf",dpi=600,pad_inches=0)
plt.show()