###############################################################################
# ReDLat plasma proteomics — machine-learning workflow
# 17. Synthesize country-held-out AUCs
# Requires: Script 12 outputs and private country metadata
# Produces: country AUC estimates, random-effects synthesis and Source Data
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

# ============================================================
# 1. IMPORTS AND USER SETTINGS
# ============================================================

import json
import math
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy import optimize, stats
from scipy.special import expit, logit
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

# ------------------------------------------------------------
# Project path
# ------------------------------------------------------------
# Leave as None to auto-detect the project folder.
PROJECT_DIR_OVERRIDE: Optional[str] = None
PROJECT_DIR = CONFIG.project_root
DATA_DIR = CONFIG.data_dir
LOCO_DIR = CONFIG.private_root / "loco" / "02_nested_cv"
OUTPUT_DIR = CONFIG.publication_root / "loco_meta_analysis"
SOURCE_DIR = OUTPUT_DIR / "source_data"
MASTER_FILE = CONFIG.master_file
OOF_FILE = LOCO_DIR / "oof_predictions.csv"
METRICS_FILE = LOCO_DIR / "metrics_nestedCV.csv"
ROC_FILE = LOCO_DIR / "roc_curves.csv"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)


def require_file(path: Path, label: str) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found:\n{path}")
    return path


for path, label in [
    (MASTER_FILE, "Master matrix"),
    (OOF_FILE, "LOCO out-of-fold predictions"),
    (METRICS_FILE, "LOCO country metrics"),
]:
    require_file(path, label)

print("PROJECT_DIR:", PROJECT_DIR)
print("MASTER_FILE:", MASTER_FILE)
print("OOF_FILE:", OOF_FILE)
print("METRICS_FILE:", METRICS_FILE)

# ============================================================
# 3. LOAD AND HARMONIZE COUNTRY-LEVEL PREDICTIONS
# ============================================================

SITE_TO_COUNTRY = {
    "BN": "Argentina",
    "BE": "Chile",
    "SL": "Chile",
    "MA": "Colombia",
    "LO": "Colombia",
    "AF": "Mexico",
    "CU": "Peru",
}

COUNTRY_ORDER = ["Argentina", "Chile", "Colombia", "Mexico", "Peru"]


def first_existing_column(
    dataframe: pd.DataFrame,
    candidates: Sequence[str],
) -> str:
    direct = {str(column): str(column) for column in dataframe.columns}
    lower = {str(column).lower(): str(column) for column in dataframe.columns}

    for candidate in candidates:
        if candidate in direct:
            return direct[candidate]
        if candidate.lower() in lower:
            return lower[candidate.lower()]

    raise KeyError(
        f"None of the expected columns were found: {list(candidates)}\n"
        f"Available columns: {list(dataframe.columns)}"
    )


master = pd.read_csv(MASTER_FILE, low_memory=False)
oof_raw = pd.read_csv(OOF_FILE)
metrics_reported = pd.read_csv(METRICS_FILE)

id_master = first_existing_column(master, ["SampleId", "sample_id", "ID"])
id_oof = first_existing_column(oof_raw, ["SampleId", "sample_id", "ID"])
target_col = first_existing_column(
    oof_raw,
    ["y_true", "target", "label", "SampleGroup_bin"],
)
score_col = first_existing_column(
    oof_raw,
    [
        "oof",
        "oof_probability",
        "oof_score",
        "decision_score",
        "y_score",
        "score",
    ],
)

master[id_master] = master[id_master].astype(str)
oof_raw[id_oof] = oof_raw[id_oof].astype(str)

metadata_columns = [id_master]
if "Country" in master.columns:
    metadata_columns.append("Country")
if "SampleGroup" in master.columns:
    metadata_columns.append("SampleGroup")

metadata = master[metadata_columns].drop_duplicates(subset=[id_master]).copy()
metadata = metadata.rename(columns={id_master: "SampleId"})

oof = oof_raw.rename(
    columns={
        id_oof: "SampleId",
        target_col: "y_true",
        score_col: "score",
    }
).copy()

oof = oof.merge(metadata, on="SampleId", how="left", validate="many_to_one")

# Country fallback from sample-ID prefix.
if "Country" not in oof.columns:
    oof["Country"] = np.nan

missing_country = oof["Country"].isna() | (oof["Country"].astype(str).str.strip() == "")
oof.loc[missing_country, "Country"] = (
    oof.loc[missing_country, "SampleId"]
    .astype(str)
    .str[:2]
    .map(SITE_TO_COUNTRY)
)

# Target fallback from SampleGroup if necessary.
oof["y_true"] = pd.to_numeric(oof["y_true"], errors="coerce")
if oof["y_true"].isna().any() and "SampleGroup" in oof.columns:
    oof.loc[oof["y_true"].isna(), "y_true"] = (
        oof.loc[oof["y_true"].isna(), "SampleGroup"]
        .map({"CN": 0, "AD": 1})
    )

oof["score"] = pd.to_numeric(oof["score"], errors="coerce")
oof["Country"] = oof["Country"].astype(str).str.strip()

analysis_data = oof.dropna(
    subset=["SampleId", "Country", "y_true", "score"]
).copy()
analysis_data["y_true"] = analysis_data["y_true"].astype(int)
analysis_data = analysis_data[
    analysis_data["Country"].isin(COUNTRY_ORDER)
].copy()

if analysis_data["SampleId"].duplicated().any():
    duplicates = analysis_data.loc[
        analysis_data["SampleId"].duplicated(keep=False),
        "SampleId",
    ].unique()
    raise ValueError(
        "Duplicate SampleId values were found in OOF predictions: "
        + ", ".join(map(str, duplicates[:10]))
    )

invalid_targets = set(analysis_data["y_true"].unique()) - {0, 1}
if invalid_targets:
    raise ValueError(f"Unexpected y_true values: {invalid_targets}")

country_counts = (
    analysis_data.groupby(["Country", "y_true"])
    .size()
    .unstack(fill_value=0)
    .rename(columns={0: "N_CN", 1: "N_AD"})
    .reset_index()
)

print("\nCountry sample sizes")
print(country_counts.to_string(index=False))
print("\nTotal participants with LOCO predictions:", len(analysis_data))

# ============================================================
# 4. DELONG AND BOOTSTRAP FUNCTIONS
# ============================================================

def compute_midrank(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    sorted_values = values[order]
    n = len(values)
    midranks = np.zeros(n, dtype=float)

    index = 0
    while index < n:
        end = index
        while end < n and sorted_values[end] == sorted_values[index]:
            end += 1
        midranks[index:end] = 0.5 * (index + end - 1) + 1
        index = end

    output = np.empty(n, dtype=float)
    output[order] = midranks
    return output


def fast_delong(
    predictions_sorted_transposed: np.ndarray,
    positive_count: int,
) -> Tuple[np.ndarray, np.ndarray]:
    m = int(positive_count)
    n = predictions_sorted_transposed.shape[1] - m
    classifiers = predictions_sorted_transposed.shape[0]

    positive_midranks = np.empty((classifiers, m))
    negative_midranks = np.empty((classifiers, n))
    all_midranks = np.empty((classifiers, m + n))

    for classifier in range(classifiers):
        positive_midranks[classifier] = compute_midrank(
            predictions_sorted_transposed[classifier, :m]
        )
        negative_midranks[classifier] = compute_midrank(
            predictions_sorted_transposed[classifier, m:]
        )
        all_midranks[classifier] = compute_midrank(
            predictions_sorted_transposed[classifier]
        )

    aucs = (
        all_midranks[:, :m].sum(axis=1) / m
        - (m + 1) / 2
    ) / n

    v01 = (
        all_midranks[:, :m] - positive_midranks
    ) / n
    v10 = 1 - (
        all_midranks[:, m:] - negative_midranks
    ) / m

    covariance_positive = np.atleast_2d(np.cov(v01, bias=False))
    covariance_negative = np.atleast_2d(np.cov(v10, bias=False))
    covariance = (
        covariance_positive / m
        + covariance_negative / n
    )
    return aucs, covariance


def delong_auc_variance(
    y_true: np.ndarray,
    scores: np.ndarray,
) -> Tuple[float, float]:
    y_true = np.asarray(y_true, dtype=int)
    scores = np.asarray(scores, dtype=float)

    if set(np.unique(y_true)) != {0, 1}:
        raise ValueError("Both CN (0) and AD (1) are required for AUC.")

    order = np.argsort(-y_true)
    positive_count = int(y_true.sum())
    predictions = scores[np.newaxis, order]

    aucs, covariance = fast_delong(
        predictions,
        positive_count,
    )
    auc_value = float(aucs[0])
    variance = float(np.atleast_2d(covariance)[0, 0])
    return auc_value, variance


def stratified_bootstrap_auc_variance(
    y_true: np.ndarray,
    scores: np.ndarray,
    repeats: int,
    random_state: int,
) -> Tuple[float, float]:
    y_true = np.asarray(y_true, dtype=int)
    scores = np.asarray(scores, dtype=float)

    positive_scores = scores[y_true == 1]
    negative_scores = scores[y_true == 0]

    if len(positive_scores) < 2 or len(negative_scores) < 2:
        raise ValueError("At least two participants per class are required.")

    rng = np.random.default_rng(random_state)
    bootstrap_aucs = np.empty(repeats, dtype=float)

    for iteration in range(repeats):
        sampled_positive = rng.choice(
            positive_scores,
            size=len(positive_scores),
            replace=True,
        )
        sampled_negative = rng.choice(
            negative_scores,
            size=len(negative_scores),
            replace=True,
        )
        bootstrap_scores = np.concatenate(
            [sampled_negative, sampled_positive]
        )
        bootstrap_labels = np.concatenate(
            [
                np.zeros(len(sampled_negative), dtype=int),
                np.ones(len(sampled_positive), dtype=int),
            ]
        )
        bootstrap_aucs[iteration] = roc_auc_score(
            bootstrap_labels,
            bootstrap_scores,
        )

    point_estimate = roc_auc_score(y_true, scores)
    variance = float(np.var(bootstrap_aucs, ddof=1))
    return float(point_estimate), variance


def safe_logit_auc(auc_value: float) -> float:
    clipped = np.clip(auc_value, 1e-6, 1 - 1e-6)
    return float(logit(clipped))


def country_auc_table(
    data: pd.DataFrame,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []

    for country_index, country in enumerate(COUNTRY_ORDER):
        country_data = data[data["Country"] == country].copy()
        if country_data.empty:
            continue

        y_true = country_data["y_true"].to_numpy(dtype=int)
        scores = country_data["score"].to_numpy(dtype=float)

        raw_auc = roc_auc_score(y_true, scores)
        score_flipped = False

        if raw_auc < 0.5 and FLIP_SCORES_IF_AUC_BELOW_HALF:
            scores = -scores
            raw_auc = roc_auc_score(y_true, scores)
            score_flipped = True

        auc_value, variance_auc = delong_auc_variance(
            y_true,
            scores,
        )
        uncertainty_method = "DeLong"

        if (
            not np.isfinite(variance_auc)
            or variance_auc <= 0
        ):
            auc_value, variance_auc = (
                stratified_bootstrap_auc_variance(
                    y_true,
                    scores,
                    repeats=BOOTSTRAP_REPEATS,
                    random_state=RANDOM_STATE + country_index,
                )
            )
            uncertainty_method = "Stratified bootstrap"

        standard_error_auc = math.sqrt(variance_auc)

        theta = safe_logit_auc(auc_value)
        derivative = 1 / (
            np.clip(auc_value, 1e-6, 1 - 1e-6)
            * np.clip(1 - auc_value, 1e-6, 1 - 1e-6)
        )
        variance_logit = variance_auc * derivative ** 2
        standard_error_logit = math.sqrt(variance_logit)

        z_critical = stats.norm.ppf(1 - ALPHA / 2)
        lower_logit = theta - z_critical * standard_error_logit
        upper_logit = theta + z_critical * standard_error_logit

        rows.append({
            "Country": country,
            "N_total": len(country_data),
            "N_CN": int((y_true == 0).sum()),
            "N_AD": int((y_true == 1).sum()),
            "AUC": auc_value,
            "SE_AUC": standard_error_auc,
            "CI95_low": float(expit(lower_logit)),
            "CI95_high": float(expit(upper_logit)),
            "Logit_AUC": theta,
            "Variance_logit_AUC": variance_logit,
            "SE_logit_AUC": standard_error_logit,
            "Uncertainty_method": uncertainty_method,
            "Score_flipped": score_flipped,
        })

    output = pd.DataFrame(rows)

    if len(output) < 2:
        raise RuntimeError(
            "At least two countries are required for synthesis."
        )

    if (output["N_CN"] == 0).any() or (output["N_AD"] == 0).any():
        raise RuntimeError(
            "Every included country must contain both CN and AD participants."
        )

    return output


country_effects = country_auc_table(analysis_data)
print("\nCountry-specific AUC estimates")
print(
    country_effects[
        [
            "Country", "N_CN", "N_AD",
            "AUC", "SE_AUC", "CI95_low", "CI95_high",
            "Uncertainty_method",
        ]
    ].to_string(index=False)
)

# ============================================================
# 5. AUDIT AGAINST THE ORIGINAL LOCO METRICS
# ============================================================

def standardize_metrics_country_column(
    dataframe: pd.DataFrame,
) -> pd.DataFrame:
    output = dataframe.copy()
    if "Country" not in output.columns:
        candidate = first_existing_column(
            output,
            ["country", "Held_out_country", "Fold"],
        )
        output = output.rename(columns={candidate: "Country"})
    return output


metrics_audit = standardize_metrics_country_column(
    metrics_reported
)

if "AUC" in metrics_audit.columns:
    metrics_audit["AUC"] = pd.to_numeric(
        metrics_audit["AUC"],
        errors="coerce",
    )
    comparison = country_effects[
        ["Country", "AUC"]
    ].merge(
        metrics_audit[["Country", "AUC"]].rename(
            columns={"AUC": "AUC_reported"}
        ),
        on="Country",
        how="left",
    )
    comparison["Absolute_difference"] = (
        comparison["AUC"] - comparison["AUC_reported"]
    ).abs()
    print("\nAUC audit against metrics_nestedCV.csv")
    print(comparison.to_string(index=False))

    if comparison["Absolute_difference"].dropna().max() > 1e-6:
        print(
            "\nWARNING: Recalculated AUCs differ from metrics_nestedCV.csv. "
            "Inspect score orientation and the OOF export."
        )
else:
    comparison = pd.DataFrame()
    print("[WARN] metrics_nestedCV.csv has no AUC column.")

# ============================================================
# 6. META-ANALYSIS FUNCTIONS
# ============================================================

def reml_tau_squared(
    effects: np.ndarray,
    variances: np.ndarray,
) -> float:
    effects = np.asarray(effects, dtype=float)
    variances = np.asarray(variances, dtype=float)

    if np.any(variances <= 0):
        raise ValueError("All sampling variances must be positive.")

    observed_variance = float(np.var(effects, ddof=1))
    upper_bound = max(
        1.0,
        observed_variance * 100,
        float(np.max(variances) * 100),
    )

    def restricted_negative_log_likelihood(
        tau_squared: float,
    ) -> float:
        total_variance = variances + tau_squared
        weights = 1 / total_variance
        pooled = np.sum(weights * effects) / np.sum(weights)
        residual = np.sum(
            weights * (effects - pooled) ** 2
        )
        value = 0.5 * (
            np.sum(np.log(total_variance))
            + np.log(np.sum(weights))
            + residual
        )
        return float(value)

    result = optimize.minimize_scalar(
        restricted_negative_log_likelihood,
        bounds=(0.0, upper_bound),
        method="bounded",
        options={"xatol": 1e-12},
    )

    if not result.success:
        raise RuntimeError(
            f"REML tau² optimization failed: {result.message}"
        )

    tau_squared = max(0.0, float(result.x))
    if tau_squared < 1e-10:
        tau_squared = 0.0
    return tau_squared


def diamond_coordinates(
    center: float,
    lower: float,
    upper: float,
    y_position: float,
    half_height: float = 0.18,
) -> np.ndarray:
    return np.array([
        [lower, y_position],
        [center, y_position + half_height],
        [upper, y_position],
        [center, y_position - half_height],
    ])


def meta_analyze_logit_auc(
    country_table: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    effects = country_table["Logit_AUC"].to_numpy(dtype=float)
    variances = country_table[
        "Variance_logit_AUC"
    ].to_numpy(dtype=float)
    k = len(effects)
    degrees_freedom = k - 1

    # Common-effect model.
    fixed_weights = 1 / variances
    fixed_effect = np.sum(fixed_weights * effects) / np.sum(fixed_weights)
    fixed_se = math.sqrt(1 / np.sum(fixed_weights))
    z_critical = stats.norm.ppf(1 - ALPHA / 2)
    fixed_lower = fixed_effect - z_critical * fixed_se
    fixed_upper = fixed_effect + z_critical * fixed_se

    # Heterogeneity based on common-effect weights.
    q_statistic = float(
        np.sum(
            fixed_weights * (effects - fixed_effect) ** 2
        )
    )
    q_p_value = float(
        stats.chi2.sf(q_statistic, degrees_freedom)
    )
    i_squared = (
        max(0.0, (q_statistic - degrees_freedom) / q_statistic)
        * 100
        if q_statistic > 0
        else 0.0
    )

    # Random-effects REML.
    tau_squared = reml_tau_squared(effects, variances)
    random_weights = 1 / (variances + tau_squared)
    random_effect = (
        np.sum(random_weights * effects)
        / np.sum(random_weights)
    )

    # Hartung-Knapp variance.
    hk_multiplier_raw = float(
        np.sum(
            random_weights
            * (effects - random_effect) ** 2
        ) / degrees_freedom
    )
    hk_multiplier = (
        max(1.0, hk_multiplier_raw)
        if MODIFIED_HARTUNG_KNAPP
        else hk_multiplier_raw
    )
    random_se_hk = math.sqrt(
        hk_multiplier / np.sum(random_weights)
    )
    t_critical = stats.t.ppf(
        1 - ALPHA / 2,
        df=degrees_freedom,
    )
    random_lower = random_effect - t_critical * random_se_hk
    random_upper = random_effect + t_critical * random_se_hk

    # Prediction interval.
    if k >= 3:
        prediction_df = k - 2
        prediction_t = stats.t.ppf(
            1 - ALPHA / 2,
            df=prediction_df,
        )
        prediction_se = math.sqrt(
            tau_squared + random_se_hk ** 2
        )
        prediction_lower = (
            random_effect - prediction_t * prediction_se
        )
        prediction_upper = (
            random_effect + prediction_t * prediction_se
        )
    else:
        prediction_lower = np.nan
        prediction_upper = np.nan

    country_output = country_table.copy()
    country_output["Weight_common_percent"] = (
        fixed_weights / fixed_weights.sum() * 100
    )
    country_output["Weight_random_percent"] = (
        random_weights / random_weights.sum() * 100
    )

    summary_rows = [
        {
            "Model": "Common effect",
            "K_countries": k,
            "Logit_AUC": fixed_effect,
            "SE_logit": fixed_se,
            "AUC": float(expit(fixed_effect)),
            "CI95_low": float(expit(fixed_lower)),
            "CI95_high": float(expit(fixed_upper)),
            "Prediction_low": np.nan,
            "Prediction_high": np.nan,
            "Tau2_logit": 0.0,
            "Tau_logit": 0.0,
            "Q": q_statistic,
            "Q_df": degrees_freedom,
            "Q_P": q_p_value,
            "I2_percent": i_squared,
            "HK_multiplier_raw": np.nan,
            "HK_multiplier_used": np.nan,
        },
        {
            "Model": "Random effects (REML + modified HK)",
            "K_countries": k,
            "Logit_AUC": random_effect,
            "SE_logit": random_se_hk,
            "AUC": float(expit(random_effect)),
            "CI95_low": float(expit(random_lower)),
            "CI95_high": float(expit(random_upper)),
            "Prediction_low": float(expit(prediction_lower)),
            "Prediction_high": float(expit(prediction_upper)),
            "Tau2_logit": tau_squared,
            "Tau_logit": math.sqrt(tau_squared),
            "Q": q_statistic,
            "Q_df": degrees_freedom,
            "Q_P": q_p_value,
            "I2_percent": i_squared,
            "HK_multiplier_raw": hk_multiplier_raw,
            "HK_multiplier_used": hk_multiplier,
        },
    ]
    summary = pd.DataFrame(summary_rows)
    return country_output, summary


country_meta, meta_summary = meta_analyze_logit_auc(
    country_effects
)

print("\nMeta-analysis summary")
print(meta_summary.to_string(index=False))

random_summary = meta_summary[
    meta_summary["Model"].str.startswith("Random")
].iloc[0]

heterogeneity_label = (
    "low"
    if random_summary["I2_percent"] < 25
    else "moderate"
    if random_summary["I2_percent"] < 50
    else "substantial"
    if random_summary["I2_percent"] < 75
    else "considerable"
)

# ============================================================
# 7. PUBLICATION-READY FOREST PLOT
# ============================================================

def format_auc_interval(
    auc_value: float,
    lower: float,
    upper: float,
) -> str:
    return f"{auc_value:.2f} [{lower:.2f}, {upper:.2f}]"


def make_forest_plot(
    country_table: pd.DataFrame,
    summary_table: pd.DataFrame,
) -> Tuple[Path, Path, Path, Path]:
    country_table = country_table.copy()
    country_table["order"] = country_table["Country"].map(
        {country: index for index, country in enumerate(COUNTRY_ORDER)}
    )
    country_table = country_table.sort_values("order").reset_index(drop=True)

    common = summary_table[
        summary_table["Model"] == "Common effect"
    ].iloc[0]
    random = summary_table[
        summary_table["Model"].str.startswith("Random")
    ].iloc[0]

    number_countries = len(country_table)
    y_country = np.arange(number_countries, 0, -1, dtype=float)
    y_common = -0.15
    y_random = -1.05
    y_prediction = -1.95

    # Compact, Nature-style panel: a dedicated label axis prevents text overlap.
    fig_width = FIGURE_WIDTH_CM / 2.54
    fig_height = FIGURE_HEIGHT_CM / 2.54
    fig = plt.figure(figsize=(fig_width, fig_height), dpi=300)
    grid = fig.add_gridspec(
        1,
        2,
        width_ratios=[1.55, 2.45],
        wspace=0.02,
    )
    label_axis = fig.add_subplot(grid[0, 0])
    axis = fig.add_subplot(grid[0, 1], sharey=label_axis)

    weights = country_table["Weight_random_percent"].to_numpy(dtype=float)
    marker_sizes = 18 + 70 * weights / weights.max()

    for row_index, row in country_table.iterrows():
        y_position = y_country[row_index]
        axis.errorbar(
            row["AUC"],
            y_position,
            xerr=np.array([
                [row["AUC"] - row["CI95_low"]],
                [row["CI95_high"] - row["AUC"]],
            ]),
            fmt="s",
            markersize=math.sqrt(marker_sizes[row_index]),
            capsize=2,
            elinewidth=0.8,
            markeredgewidth=0.6,
        )
        axis.text(
            1.012,
            y_position,
            f"{row['AUC']:.2f}",
            fontsize=5.7,
            ha="left",
            va="center",
            clip_on=False,
        )

    common_diamond = Polygon(
        diamond_coordinates(
            common["AUC"],
            common["CI95_low"],
            common["CI95_high"],
            y_common,
        ),
        closed=True,
    )
    axis.add_patch(common_diamond)

    random_diamond = Polygon(
        diamond_coordinates(
            random["AUC"],
            random["CI95_low"],
            random["CI95_high"],
            y_random,
        ),
        closed=True,
    )
    axis.add_patch(random_diamond)

    axis.errorbar(
        random["AUC"],
        y_prediction,
        xerr=np.array([
            [random["AUC"] - random["Prediction_low"]],
            [random["Prediction_high"] - random["AUC"]],
        ]),
        fmt="o",
        markersize=3.2,
        capsize=2,
        elinewidth=0.8,
    )

    axis.axvline(
        0.5,
        linestyle="--",
        linewidth=0.8,
        alpha=0.7,
    )
    axis.set_xlim(*X_LIMITS)
    axis.set_ylim(y_prediction - 0.85, number_countries + 0.85)
    axis.set_yticks([])
    axis.set_xlabel("ROC AUC", fontsize=7)
    axis.set_xticks([0.50, 0.60, 0.70, 0.80, 0.90, 1.00])
    axis.tick_params(axis="x", labelsize=6, length=2, width=0.5)
    axis.set_title(
        "Country-held-out ROC AUC synthesis",
        fontsize=7.4,
        pad=7,
    )

    # Left text axis.
    label_axis.set_xlim(0, 1)
    label_axis.set_ylim(y_prediction - 0.85, number_countries + 0.85)
    label_axis.set_xticks([])
    label_axis.set_yticks([])
    for spine in label_axis.spines.values():
        spine.set_visible(False)

    label_axis.text(
        0.00,
        number_countries + 0.52,
        "Held-out country",
        fontsize=5.7,
        fontweight="bold",
        ha="left",
        va="center",
    )
    label_axis.text(
        0.98,
        number_countries + 0.52,
        "CN / AD",
        fontsize=5.7,
        fontweight="bold",
        ha="right",
        va="center",
    )

    for row_index, row in country_table.iterrows():
        y_position = y_country[row_index]
        label_axis.text(
            0.00,
            y_position,
            row["Country"],
            fontsize=6,
            ha="left",
            va="center",
        )
        label_axis.text(
            0.98,
            y_position,
            f"{int(row['N_CN'])} / {int(row['N_AD'])}",
            fontsize=5.7,
            ha="right",
            va="center",
        )

    for label, y_position in [
        ("Common effect", y_common),
        ("Random effects", y_random),
        ("95% prediction interval", y_prediction),
    ]:
        label_axis.text(
            0.00,
            y_position,
            label,
            fontsize=5.8 if "prediction" in label else 6,
            fontweight="normal" if "prediction" in label else "bold",
            ha="left",
            va="center",
        )

    axis.text(
        1.012,
        y_common,
        f"{common['AUC']:.2f}",
        fontsize=5.8,
        fontweight="bold",
        ha="left",
        va="center",
        clip_on=False,
    )
    axis.text(
        1.012,
        y_random,
        f"{random['AUC']:.2f}",
        fontsize=5.8,
        fontweight="bold",
        ha="left",
        va="center",
        clip_on=False,
    )

    heterogeneity_text = (
        f"Q = {random['Q']:.2f}, df = {int(random['Q_df'])}, "
        f"P = {random['Q_P']:.3g}; I² = {random['I2_percent']:.1f}%; "
        f"tau² = {random['Tau2_logit']:.3f}"
    )
    fig.text(
        0.04,
        0.055,
        heterogeneity_text,
        fontsize=5.4,
        ha="left",
        va="center",
    )

    label_axis.text(
        -0.18,
        1.08,
        PANEL_LETTER,
        transform=label_axis.transAxes,
        fontsize=10,
        fontweight="bold",
        ha="left",
        va="top",
    )

    axis.spines["left"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["top"].set_visible(False)
    axis.spines["bottom"].set_linewidth(0.6)

    fig.subplots_adjust(
        left=0.04,
        right=0.90,
        top=0.86,
        bottom=0.19,
    )

    pdf_path = OUTPUT_DIR / "Figure4e_LOCO_AUC_meta_analysis.pdf"
    svg_path = OUTPUT_DIR / "Figure4e_LOCO_AUC_meta_analysis.svg"
    png_path = OUTPUT_DIR / "Figure4e_LOCO_AUC_meta_analysis.png"
    tiff_path = OUTPUT_DIR / "Figure4e_LOCO_AUC_meta_analysis.tiff"

    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(png_path, bbox_inches="tight", dpi=600)
    fig.savefig(
        tiff_path,
        bbox_inches="tight",
        dpi=600,
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.show()

    return pdf_path, svg_path, png_path, tiff_path


figure_paths = make_forest_plot(
    country_meta,
    meta_summary,
)

print("\nFigure files")
for path in figure_paths:
    print(path)

# ============================================================
# 8. SOURCE DATA AND MACHINE-READABLE AUDIT
# ============================================================

country_csv = SOURCE_DIR / "LOCO_country_AUC_estimates.csv"
summary_csv = SOURCE_DIR / "LOCO_AUC_meta_analysis_summary.csv"
audit_csv = SOURCE_DIR / "LOCO_AUC_recalculation_audit.csv"
source_xlsx = SOURCE_DIR / "LOCO_AUC_meta_analysis_source_data.xlsx"
manifest_json = SOURCE_DIR / "LOCO_AUC_meta_analysis_manifest.json"

country_meta.to_csv(country_csv, index=False)
meta_summary.to_csv(summary_csv, index=False)
if not comparison.empty:
    comparison.to_csv(audit_csv, index=False)

with pd.ExcelWriter(source_xlsx, engine="openpyxl") as writer:
    country_meta.to_excel(
        writer,
        sheet_name="Country_AUC",
        index=False,
    )
    meta_summary.to_excel(
        writer,
        sheet_name="Meta_summary",
        index=False,
    )
    country_counts.to_excel(
        writer,
        sheet_name="Country_counts",
        index=False,
    )
    if not comparison.empty:
        comparison.to_excel(
            writer,
            sheet_name="AUC_audit",
            index=False,
        )

manifest = {
    "analysis_name": "Country-held-out ROC AUC synthesis",
    "interpretation": (
        "Internal geographic robustness; not external validation."
    ),
    "project_dir": str(PROJECT_DIR),
    "inputs": {
        "master_file": str(MASTER_FILE),
        "oof_predictions": str(OOF_FILE),
        "reported_metrics": str(METRICS_FILE),
        "roc_curves": str(ROC_FILE),
    },
    "settings": {
        "alpha": ALPHA,
        "uncertainty": (
            "DeLong variance with stratified-bootstrap fallback"
        ),
        "effect_scale": "logit AUC",
        "tau2_estimator": "REML",
        "inference": (
            "modified Hartung-Knapp"
            if MODIFIED_HARTUNG_KNAPP
            else "Hartung-Knapp"
        ),
        "bootstrap_repeats": BOOTSTRAP_REPEATS,
        "random_state": RANDOM_STATE,
        "flip_scores_below_half": FLIP_SCORES_IF_AUC_BELOW_HALF,
    },
    "outputs": {
        "country_csv": str(country_csv),
        "summary_csv": str(summary_csv),
        "source_xlsx": str(source_xlsx),
        "figure_pdf": str(figure_paths[0]),
        "figure_svg": str(figure_paths[1]),
        "figure_png": str(figure_paths[2]),
        "figure_tiff": str(figure_paths[3]),
    },
}

with open(
    manifest_json,
    "w",
    encoding="utf-8",
) as file:
    json.dump(
        manifest,
        file,
        indent=2,
        ensure_ascii=False,
    )

print("\nSource data")
print(country_csv)
print(summary_csv)
print(source_xlsx)
print(manifest_json)

# ============================================================
# 9. GENERATE METHODS, RESULTS AND LEGEND TEXT
# ============================================================

minimum_row = country_meta.loc[country_meta["AUC"].idxmin()]
maximum_row = country_meta.loc[country_meta["AUC"].idxmax()]

methods_text = (
    "Inter-country synthesis of classification performance. "
    "Country-specific discrimination was estimated from predictions generated "
    "exclusively in the held-out country during leave-one-country-out evaluation. "
    "ROC AUCs and DeLong standard errors were calculated for each country from "
    "the corresponding held-out predictions. AUCs were logit-transformed and "
    "synthesized using inverse-variance common-effect and random-effects models. "
    "Between-country variance was estimated by restricted maximum likelihood, "
    "and random-effects inference used a modified Hartung-Knapp adjustment. "
    "Summary estimates and confidence and prediction intervals were "
    "back-transformed to the AUC scale. Heterogeneity was quantified using "
    "Cochran's Q, I² and tau². Because all countries belonged to ReDLat and LOCO "
    "training sets partially overlapped, this analysis was interpreted as a "
    "descriptive synthesis of internal geographic robustness rather than "
    "external validation."
)

results_text = (
    f"Country-held-out AUCs ranged from {minimum_row['AUC']:.2f} in "
    f"{minimum_row['Country']} to {maximum_row['AUC']:.2f} in "
    f"{maximum_row['Country']}. Random-effects synthesis yielded a pooled AUC "
    f"of {random_summary['AUC']:.2f} "
    f"(95% CI, {random_summary['CI95_low']:.2f}-"
    f"{random_summary['CI95_high']:.2f}; 95% prediction interval, "
    f"{random_summary['Prediction_low']:.2f}-"
    f"{random_summary['Prediction_high']:.2f}), with {heterogeneity_label} "
    f"between-country heterogeneity (I² = "
    f"{random_summary['I2_percent']:.1f}%). These findings indicate that "
    "discrimination was retained across held-out recruitment settings, while "
    "its magnitude varied geographically."
)

legend_text = (
    "e, Forest plot showing country-specific ROC AUCs obtained during "
    "leave-one-country-out evaluation. Squares represent country-level AUCs, "
    "horizontal lines indicate 95% confidence intervals and square size is "
    "proportional to random-effects inverse-variance weight. Diamonds show the "
    "common-effect and random-effects pooled estimates. The prediction interval "
    "indicates the expected range of performance in a comparable recruitment "
    "setting. AUCs were synthesized on the logit scale using restricted maximum "
    "likelihood and modified Hartung-Knapp inference. This analysis assesses "
    "internal geographic robustness and does not constitute external validation."
)

report_path = SOURCE_DIR / "LOCO_meta_analysis_reporting_text.txt"
with open(report_path, "w", encoding="utf-8") as file:
    file.write("METHODS\n")
    file.write(methods_text + "\n\n")
    file.write("RESULTS\n")
    file.write(results_text + "\n\n")
    file.write("FIGURE LEGEND\n")
    file.write(legend_text + "\n")

print("\nMETHODS\n")
print(methods_text)
print("\nRESULTS\n")
print(results_text)
print("\nFIGURE LEGEND\n")
print(legend_text)
print("\nSaved reporting text:", report_path)

# ============================================================
# 10. FINAL SCIENTIFIC CHECKS
# ============================================================

checks = pd.DataFrame([
    {
        "Check": "Five countries included",
        "Pass": len(country_meta) == 5,
        "Observed": len(country_meta),
    },
    {
        "Check": "Every country has CN participants",
        "Pass": bool((country_meta["N_CN"] > 0).all()),
        "Observed": int(country_meta["N_CN"].min()),
    },
    {
        "Check": "Every country has AD participants",
        "Pass": bool((country_meta["N_AD"] > 0).all()),
        "Observed": int(country_meta["N_AD"].min()),
    },
    {
        "Check": "All AUCs at or above chance",
        "Pass": bool((country_meta["AUC"] >= 0.5).all()),
        "Observed": float(country_meta["AUC"].min()),
    },
    {
        "Check": "Random-effects CI within valid AUC scale",
        "Pass": bool(
            0 <= random_summary["CI95_low"]
            <= random_summary["CI95_high"]
            <= 1
        ),
        "Observed": (
            f"{random_summary['CI95_low']:.3f}-"
            f"{random_summary['CI95_high']:.3f}"
        ),
    },
])

checks_path = SOURCE_DIR / "LOCO_meta_analysis_final_checks.csv"
checks.to_csv(checks_path, index=False)

print("\nFinal checks")
print(checks.to_string(index=False))

if not checks["Pass"].all():
    raise RuntimeError(
        "At least one final scientific check failed. "
        f"Inspect: {checks_path}"
    )

print("\nAnalysis completed successfully.")


