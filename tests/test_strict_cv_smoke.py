from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from redlat_ml.strict_cv import (
    StrictCVSettings, select_stable_panel, fit_tuned_calibrated_svm,
    predict_probability, metric_row, fold_permutation_importance,
)

rng = np.random.default_rng(42)
n = 90
p = 24
y = pd.Series(np.r_[np.zeros(n//2, dtype=int), np.ones(n//2, dtype=int)])
X = pd.DataFrame(rng.normal(size=(n, p)), columns=[f'P{i:02d}' for i in range(p)])
X.loc[y.eq(1), 'P00'] += 1.2
X.loc[y.eq(1), 'P01'] -= 0.8
missing = rng.random(X.shape) < 0.05
X = X.mask(missing)
train = np.r_[np.arange(0, 36), np.arange(45, 81)]
test = np.r_[np.arange(36, 45), np.arange(81, 90)]
settings = StrictCVSettings(
    random_state=42, inner_folds=3, stability_threshold=0.66,
    selector_search_iter=2, svm_search_iter=2, permutation_repeats=2,
    n_jobs=1, fallback_top_n=5,
)
selection = select_stable_panel(X.iloc[train], y.iloc[train], X.columns, settings, 42)
assert selection['panel']
estimator, params = fit_tuned_calibrated_svm(X.iloc[train][selection['panel']], y.iloc[train], settings, 42)
prob = predict_probability(estimator, X.iloc[test][selection['panel']])
assert np.isfinite(prob).all() and ((prob >= 0) & (prob <= 1)).all()
row = metric_row(y.iloc[test], prob)
assert np.isfinite(row['AUC'])
importance = fold_permutation_importance(estimator, X.iloc[test][selection['panel']], y.iloc[test], selection['panel'], settings, 'smoke')
assert not importance.empty
print({'panel_size': len(selection['panel']), 'auc': row['AUC'], 'params': params})
