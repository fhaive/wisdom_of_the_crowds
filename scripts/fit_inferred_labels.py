import numpy as np
import pandas as pd
import xarray as xa
from functools import reduce
from wisdom.data import load_physchem_data
from wisdom.data import load_gene_expression
from wisdom.training import generate_splits
from wisdom.training import early_multiview_crossvalidate
from wisdom.training import late_multiview_crossvalidate
from wisdom.training import crossval_evaluate
from scipy.stats import gaussian_kde
from loguru import logger
from scipy.stats import beta
from scipy.stats import gaussian_kde
import click

@click.command()
@click.option("--physchem_properties", type=click.Path(dir_okay=False, exists=True), default="data/physicochemical_descriptors.txt")
@click.option("--gex_phenodata", type=click.Path(dir_okay=False, exists=True), default="data/phenodata.txt")
@click.option("--gex", type=click.Path(dir_okay=False, exists=True), default="data/gex.csv.gz")
@click.option("--inferred_labels", type=click.Path(dir_okay=False, exists=True), default="outputs/model_inference_anon.nc")
@click.option("--n_repeats", type=int, default=100)
@click.option("--n_folds", type=int, default=5)
@click.option("--topk_features", type=int, default=10)
@click.option("--n_processes", type=int, default=32)
@click.option("--output_folder", type=click.Path(file_okay=False, exists=True), default="outputs/cross_validation")
def main(physchem_properties, gex_phenodata, gex, inferred_labels, output_folder, n_repeats, n_folds, topk_features, n_processes):
    INFERENCE_FILE = inferred_labels
    PHYSCHEM_FILE = physchem_properties
    GEX_PHENODATA_FILE = gex_phenodata
    GEX_DATA_FILE = gex
    OUTPUT_FOLDER = output_folder

    N_REPEATS = n_repeats
    N_FOLDS = n_folds
    TOPK_FEATURES = topk_features
    N_PROCESSES = n_processes

    ENDPOINTS = [
        "Asthma",
        "Carcinogenic",
        "Cardiovascular.toxicity",
        "Cytotoxic",
        "Deregulation.of.cellular.metabolism",
        "Endocrine.deregulation.disruption",
        "Environmental.hazard",
        "Genotoxic",
        "Hepatotoxic",
        "Immunotoxic",
        "Intestinal.toxicity",
        "Lung.fibrosis",
        "Nephrotoxic",
        "Neurotoxic",
        "Reproductive.toxicity",
        "Skin.irritation",
        "Skin.sensitizatioin",
        "Global",
    ]

    logger.info("Loading estimated labels")
    concerns = xa.open_dataset(INFERENCE_FILE)
    labels, weights = labels_and_uncertainties(concerns)
    print(weights.shape)


    logger.info("Loading physicochemical properties")
    X_physchem = load_physchem_data(PHYSCHEM_FILE)
    logger.info("Loading gene expression aggregated by median")
    X_gex_median = load_gene_expression(
        GEX_PHENODATA_FILE,
        GEX_DATA_FILE,
        aggregate="median",
        min_std_thresh=0,
    )
    logger.info("Loading gene expression aggregated by absolute logFC")
    X_gex_absmax = load_gene_expression(
        GEX_PHENODATA_FILE,
        GEX_DATA_FILE,
        aggregate="abs_max",
        min_std_thresh=0
    )

    X_gex_absmax, X_gex_median, X_physchem, labels, weights = filter_common_elements(
        X_gex_absmax, X_gex_median, X_physchem, labels, weights
    )

    logger.debug("Labels shape {}".format(labels.shape))
    logger.debug("Weights shape {}".format(weights.shape))
    logger.debug("Physchem shape {}".format(X_physchem.shape))
    logger.debug("Gex median shape {}".format(X_gex_median.shape))
    logger.debug("Gex absmax shape {}".format(X_gex_absmax.shape))

    for e in ENDPOINTS:
        logger.info("Analysing {} endpoint".format(e))
        y = labels[e]
        w = weights[e]
        common_splits = list(generate_splits(y, n_repeats=N_REPEATS, n_folds=N_FOLDS))

        logger.info("Fitting single view models on gex_absmax data")
        gex_absmax_perf_df = crossval_evaluate(
            X_gex_absmax,
            y,
            w,
            n_repeats=N_REPEATS,
            n_folds=N_FOLDS,
            n_processes=N_PROCESSES,
            splits=common_splits
        )
        gex_absmax_perf_df.to_json(OUTPUT_FOLDER + "gex_absmax_{}.json".format(e))

        logger.info("Fitting single view models on gex_median data")
        gex_median_perf_df = crossval_evaluate(
            X_gex_median,
            y,
            w,
            n_repeats=N_REPEATS,
            n_folds=N_FOLDS,
            n_processes=N_PROCESSES,
            splits=common_splits
        )
        gex_median_perf_df.to_json(OUTPUT_FOLDER + "gex_median_{}.json".format(e))

        logger.info("Fitting single view models on physchem data")
        physchem_perf_df = crossval_evaluate(
            X_physchem,
            y,
            w,
            n_repeats=N_REPEATS,
            n_folds=N_FOLDS,
            n_processes=N_PROCESSES,
            splits=common_splits
        )
        physchem_perf_df.to_json(OUTPUT_FOLDER + "physchem_{}.json".format(e))

        logger.info("Fitting top 10 gex_median + physchem early multi view model")
        early_top10_median_multi_view_perf_df = early_multiview_crossvalidate(
            X_physchem,
            X_gex_median,
            y,
            w,
            physchem_perf_df,
            gex_median_perf_df,
            common_splits,
            top=TOPK_FEATURES,
            n_processes=N_PROCESSES
        )
        early_top10_median_multi_view_perf_df.to_json(OUTPUT_FOLDER + "early-top10_median_{}.json".format(e))

        logger.info("Fitting top 10 gex_absmax + physchem early multi view model")
        early_top10_absmax_multi_view_perf_df = early_multiview_crossvalidate(
            X_physchem,
            X_gex_absmax,
            y,
            w,
            physchem_perf_df,
            gex_median_perf_df,
            common_splits,
            top=TOPK_FEATURES,
            n_processes=N_PROCESSES
        )
        early_top10_absmax_multi_view_perf_df.to_json(OUTPUT_FOLDER + "early-top10_absmax_{}.json".format(e))

        logger.info("Fitting gex_median + physchem early multi view model")
        early_median_multi_view_perf_df = early_multiview_crossvalidate(
            X_physchem,
            X_gex_median,
            y,
            w,
            physchem_perf_df,
            gex_median_perf_df,
            common_splits,
            top=None,
            n_processes=N_PROCESSES
        )
        early_median_multi_view_perf_df.to_json(OUTPUT_FOLDER + "early_median_{}.json".format(e))

        logger.info("Fitting gex_absmax + physchem early multi view model")
        early_absmax_multi_view_perf_df = early_multiview_crossvalidate(
            X_physchem,
            X_gex_absmax,
            y,
            w,
            physchem_perf_df,
            gex_median_perf_df,
            common_splits,
            top=None,
            n_processes=N_PROCESSES
        )
        early_absmax_multi_view_perf_df.to_json(OUTPUT_FOLDER + "early_absmax_{}.json".format(e))

        logger.info("Fitting late multi view models")
        late_multi_view_perf_df = late_multiview_crossvalidate(physchem_perf_df, gex_median_perf_df, w)
        late_multi_view_perf_df.to_json(OUTPUT_FOLDER + "late_{}.json".format(e))

def mode(data, num_samples=1000):
    if len(data) == 1:
        return data.values[0]
    xmin = np.min(data)
    xmax = np.max(data)
    x = np.linspace(xmin, xmax, num_samples)
    y = gaussian_kde(data)(x)
    return x[y.argmax()]

def binary_entropy(a, axis=0):
    n = a.shape[axis]
    p_one = a.sum(axis=axis) / n
    p_zero = 1 - p_one
    return -(p_zero * np.log2(p_zero) + p_one * np.log2(p_one))

def labels_and_uncertainties(inference):
    y = inference.inferred_labels.to_pandas().T
    # compute sample weights for all the hazard endpoints
    samples = inference.samples.to_numpy()
    H = binary_entropy(samples, axis=0).T
    # values closer to 1 (less entropy) should be taken into account more by the classifiers
    w = pd.DataFrame(1 - H, index=y.index, columns=y.columns)
    # use estimated parameters of the posterior concerns to estimate weights for the global concern
    p_a = inference.posterior_concerns_a.to_pandas()
    p_b = inference.posterior_concerns_b.to_pandas()
    y["Global"] = (inference.posterior_global_concerns.to_pandas().apply(mode) > 0.5).astype(float)
    # since differential entropy is negative we transform it using 1-exp(H) to be in [0, 1], where larger values are associated to more certain concerns
    w["Global"] = (1-np.exp(pd.Series(beta(p_a, p_b).entropy(), index=p_a.index)))
    return y, w

def filter_common_elements(*args):
    common_elements = reduce(
        set.intersection, map(lambda x: set(x.index), args)
    )
    return tuple(x.reindex(common_elements) for x in args)

if __name__ == "__main__":
    main()