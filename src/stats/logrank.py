import os
import sys
import numpy as np
import pandas as pd
import lifelines

from tqdm import tqdm

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath("__file__"))))
from src.stats.rmst import calculate_RMST_for_one_dataset


def calculate_statistical_test(
    tcga_datasets: list,
    path_to_clinical: str = "data/raw",
    geneset_name: str = "Hallmarks",
    return_value: str = "p-value",
    rmst_time_limit: float = 0.75,
) -> pd.DataFrame:
    results = []

    for tcga in tqdm(tcga_datasets):
        metadata = pd.read_csv(
            f"{path_to_clinical}/TCGA-{tcga}-survival.csv",
            sep=",",
            index_col=["samples"],
        )
        ssgsea = pd.read_csv(
            f"data/{geneset_name}/TCGA-{tcga}-ssGSEA.tsv",
            sep="\t",
            index_col=["samples"],
        )
        ssgsea = ssgsea.loc[ssgsea.index.drop_duplicates(keep=False)]

        common_ids = list(set(metadata.index.values).intersection(ssgsea.index.values))
        metadata = metadata.loc[common_ids]
        ssgsea = ssgsea.loc[common_ids]

        if return_value != "rmst":
            results += [
                calculate_logrank_for_one_dataset(metadata, ssgsea, tcga, return_value)
            ]
        else:
            results += [
                calculate_RMST_for_one_dataset(
                    metadata, ssgsea, tcga, rmst_time_limit=rmst_time_limit
                )
            ]

    return pd.concat(results, axis=1)


def calculate_logrank_for_one_dataset(
    metadata: pd.DataFrame,
    ssgsea: pd.DataFrame,
    tcga: str,
    return_value: str = "p-value",
) -> pd.Series:
    scores = {}

    for split_by in ssgsea.columns:
        up_ids = ssgsea[split_by] > np.median(ssgsea[split_by])
        up_score = metadata.loc[up_ids[up_ids].index]
        down_score = metadata.loc[up_ids[~up_ids].index]

        logrank_test_score = lifelines.statistics.logrank_test(
            up_score["time"],
            down_score["time"],
            event_observed_A=up_score["event"],
            event_observed_B=down_score["event"],
        )

        if return_value == "p-value":
            scores[split_by] = logrank_test_score.p_value
        elif return_value == "test-statistic":
            scores[split_by] = logrank_test_score.test_statistic
        else:
            raise KeyError(
                f"Invalid return_value: {return_value}. Choose either p-value or test-statistic."
            )

    return pd.Series(scores, name=tcga)
