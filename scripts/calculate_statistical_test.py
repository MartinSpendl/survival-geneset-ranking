import os
import sys

sys.path.append(os.path.dirname(os.path.abspath("__file__")))
from src.utils.constants import TCGA_DATASETS
from src.stats.logrank import calculate_statistical_test

if __name__ == "__main__":

    PATH_TO_CLINICAL = "tcga-data/data/raw"
    GENESET_NAME = "Hallmarks"  # "GO-BiologicalProcess"
    RMST_TIME_LIMIT = 0.75

    # Logrank p-value
    logrank_score = calculate_statistical_test(
        TCGA_DATASETS, PATH_TO_CLINICAL, GENESET_NAME, return_value="p-value"
    )

    logrank_score.reset_index().to_csv(
        f"data/{GENESET_NAME}/logrank-scores-p-value.tsv", sep="\t", index=False
    )

    # Logrank test-statistic
    logrank_score = calculate_statistical_test(
        TCGA_DATASETS, PATH_TO_CLINICAL, GENESET_NAME, return_value="test-statistic"
    )

    logrank_score.reset_index().to_csv(
        f"data/{GENESET_NAME}/logrank-scores-test-statistic.tsv", sep="\t", index=False
    )

    # RMST difference
    rmst_score = calculate_statistical_test(
        TCGA_DATASETS,
        PATH_TO_CLINICAL,
        GENESET_NAME,
        return_value="rmst",
        rmst_time_limit=RMST_TIME_LIMIT,
    )

    rmst_score.reset_index().to_csv(
        f"data/{GENESET_NAME}/rmst-difference-scores.tsv", sep="\t", index=False
    )
