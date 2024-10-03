import numpy as np
import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
from xgboost.sklearn import XGBClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score

def crossval_evaluate(X, y, weights=None, n_repeats=10, n_folds=5, n_processes=6, splits=None):
    tqdm._instances.clear()
    with mp.Pool(processes=n_processes) as p:
        if splits is None:
            splits = generate_splits(y, n_repeats, n_folds)
        splits_with_data = (s + (X, y, weights) for s in splits)
        perf = list(tqdm(p.imap_unordered(xgboost_fit_and_score, splits_with_data)))
    perf_df = pd.DataFrame(
        perf,
        columns=[
            "repeat",
            "fold",
            "train_idx",
            "test_idx",
            "train_score",
            "test_score",
            "train_pred",
            "test_pred",
            "test_true",
            "feature_importances",
        ],
    )
    return perf_df

def late_multiview_crossvalidate(physchem_perf, gex_perf, weights=None):
    repeats = sorted(list(physchem_perf["repeat"].unique()))
    folds = sorted(list(physchem_perf["fold"].unique()))

    physchem_repeats = aggregate_test_predictions(physchem_perf).groupby("repeat")
    gex_repeats = aggregate_test_predictions(gex_perf).groupby("repeat")

    stacked_perf_df = []

    for r in repeats:
        true_labels = physchem_repeats.get_group(r).set_index("test_name")["test_true"]
        physchem_preds = physchem_repeats.get_group(r).set_index("test_name")["test_pred"]
        gex_preds = gex_repeats.get_group(r).set_index("test_name")["test_pred"]
        X = pd.concat(
            [
                pd.DataFrame(
                    gex_preds.to_list(), index=gex_preds.index
                ).add_suffix("_gex"),
                pd.DataFrame(
                    physchem_preds.to_list(), index=physchem_preds.index
                ).add_suffix("_physchem"),
                true_labels,
            ],
            axis=1,
            join="inner",
        )

        y = X.pop("test_true").astype(float)

        perf_df = crossval_evaluate(X, y, weights, n_repeats=1, n_folds=5)
        perf_df["repeat"] = r
        stacked_perf_df.append(perf_df)
    stacked_perf_df = pd.concat(stacked_perf_df, axis=0, ignore_index=True)
    return stacked_perf_df

def early_multiview_crossvalidate(
    X_physchem, X_gex, y, weights, physchem_perf, gex_perf, common_splits, top=None, n_processes=4
):
    tqdm._instances.clear()
    X = pd.merge(X_physchem, X_gex, left_index=True, right_index=True, suffixes=("_gex", "_physchem"))
    with mp.Pool(processes=n_processes) as p:
        if top is not None:
            splits_with_data = (
                s 
                + (
                    X[get_top_features(
                        X_physchem.columns,
                        X_gex.columns,
                        physchem_perf,
                        gex_perf,
                        top,
                        s
                    )],
                    y,
                    weights
                ) for s in common_splits
            )
        else:
            splits_with_data = (s + (X, y, weights) for s in common_splits)
        perf = list(tqdm(p.imap_unordered(xgboost_fit_and_score, splits_with_data)))
    perf_df = pd.DataFrame(
        perf,
        columns=[
            "repeat",
            "fold",
            "train_idx",
            "test_idx",
            "train_score",
            "test_score",
            "train_pred",
            "test_pred",
            "test_true",
            "feature_importances",
        ],
    )
    return perf_df

def get_top_features(physchem_features, gex_features, physchem_perf, gex_perf, top, split):
    top_physchem = physchem_features[
        np.argpartition(
            physchem_perf.loc[
                (physchem_perf["repeat"]==split[0])
                & (physchem_perf["fold"]==split[1]), 
                "feature_importances"
            ].values[0],
            -top
        )[-top:]
    ].to_list()
    top_gex = gex_features[
        np.argpartition(
            gex_perf.loc[
                (gex_perf["repeat"]==split[0])
                & (gex_perf["fold"]==split[1]), 
                "feature_importances"
            ].values[0],
            -top
        )[-top:]
    ].to_list()
    return top_physchem + top_gex

def generate_splits(y, n_repeats=10, n_folds=5):
    names = np.array(y.index.to_list())
    for repeat in range(n_repeats):
        for i, (train_idx, test_idx) in enumerate(
            StratifiedKFold(n_splits=n_folds, shuffle=True).split(y, y)
        ):
            train_names = [names[i] for i in train_idx]
            test_names = [names[i] for i in test_idx]
            # X_train, X_test = X[train_idx], X[test_idx]
            # y_train, y_test = y[train_idx], y[test_idx]
            # yield (X_train, X_test, y_train, y_test)
            yield (repeat, i, train_names, test_names)

def xgboost_fit_and_score(args):
    repeat, fold, train_idx, test_idx, X, y, weights = args
    X_train, X_test = X.loc[train_idx], X.loc[test_idx]
    y_train, y_test = y.loc[train_idx], y.loc[test_idx]
    if weights is not None:
        w_train = weights.loc[train_idx]
    else:
        w_train = None
    rf = XGBClassifier(
        objective="binary:logistic",
        n_estimators=5000,
        importance_type="gain",
        eval_metric="logloss",
        max_depth=50,
        n_jobs=1,
        missing=np.nan,
    )
    rf.fit(X_train.values, y_train.values, sample_weight=w_train)
    y_train_pred = rf.predict_proba(X_train.values)
    y_test_pred = rf.predict_proba(X_test.values)
    return (
        repeat,
        fold,
        train_idx,
        test_idx,
        accuracy_score(y_train, y_train_pred.argmax(axis=1)),
        accuracy_score(y_test, y_test_pred.argmax(axis=1)),
        # rf.oob_score_,
        y_train_pred,
        y_test_pred,
        y_test,
        rf.feature_importances_,
    )


def aggregate_test_predictions(perf_df):
    test_name = perf_df["test_true"].apply(
        lambda x: x.index.to_list()
    )
    preds = perf_df.loc[:, ["repeat", "fold", "test_pred", "test_true"]]
    preds["test_name"] = test_name
    preds = preds.explode(column=["test_pred", "test_true", "test_name"])
    return preds
