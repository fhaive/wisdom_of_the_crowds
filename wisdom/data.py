import numpy as np
import pandas as pd

def stack_table(table):
    df = (
        table.stack()
        .reset_index()
        .rename(columns={"level_1": "hazard", 0: "prediction"})
        )
    return df


def load_responses(datadir):
    crowd = pd.read_excel(
        datadir + "Combined_cleaned_responses_anon.xlsx",
        sheet_name=None,
        na_values="x",
        false_values="n",
        true_values="y",
        index_col=0,
    )
    table_list = []
    for key, table in crowd.items():
        df = stack_table(table)
        df["expert"] = key
        table_list.append(df)

    data = pd.concat(table_list).pivot(
        index=["ENM", "hazard"], columns="expert", values="prediction"
    ).astype(np.float32)
    return data

def load_gene_expression(
    pheno_file,
    data_file,
    aggregate="median",
    max_nan_ratio=0.1,
    min_std_thresh=0.5
):
    pheno = pd.read_csv(pheno_file, index_col=0)
    pheno = pheno.loc[pheno["in.vivo-in.vitro"]=="in vitro"]
    gex = pd.read_csv(data_file, index_col=0)
    nan_ratio = gex.isna().sum(axis=0) / gex.shape[0]
    gex = gex.loc[:, nan_ratio <= max_nan_ratio]
    gex_group = gex.groupby(pheno["woc_sample"])

    if aggregate == "median":
        gex = gex_group.median()
    elif aggregate == "abs_max":
        df_max = gex_group.max()
        df_min = gex_group.min()
        gex = pd.DataFrame(
            np.where(df_max > -df_min, df_max, df_min),
            index=df_max.index,
            columns=df_max.columns
        )
    else:
        ValueError("Unknown aggregate method: {}".format(aggregate))

    gex = gex.loc[:, gex.std() > min_std_thresh]
    return gex

def load_physchem_data(data_file):
    PHYSCHEM_TO_DROP = [
        'Nominal_diameter',
        'Nominal_length',
        'Nominal_Specific_Surface_area',
        'TEM_diameter',
        'TEM_width',
        'TEM_length',
        'BET_surface_area',
        'DLS_Mean_Diameter_water',
        'DLS_Mean_Diameter_medium',
        'Zeta_Potential_water',
        'Zeta_Potential_medium',
        "purity",
        "Endotoxins",
    ]

    physchem = pd.read_csv(data_file, index_col=0)
    return physchem.drop(columns=PHYSCHEM_TO_DROP)
