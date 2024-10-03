import numpy as np
from scipy.stats import gaussian_kde

def aggregate_test_predictions(perf_df):
    test_name = perf_df["test_true"].apply(
        lambda x: x.index.to_list()
    )
    preds = perf_df.loc[:, ["repeat", "fold", "test_pred", "test_true"]]
    preds["test_name"] = test_name
    preds = preds.explode(column=["test_pred", "test_true", "test_name"])
    return preds

def mode(data, num_samples=1000):
    xmin = np.min(data)
    xmax = np.max(data)
    x = np.linspace(xmin, xmax, num_samples)
    y = gaussian_kde(data)(x)
    return x[y.argmax()]

def empyrical_mode(x, bins=20):
    counts, bins = np.histogram(x, bins=bins)
    max_bin = np.argmax(counts)
    return bins[max_bin:max_bin+2].mean()