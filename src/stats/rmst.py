import numpy as np
import pandas as pd


def calculate_RMST_for_one_dataset(
    metadata: pd.DataFrame,
    ssgsea: pd.DataFrame,
    tcga: str,
    rmst_time_limit: float = 0.75,
) -> pd.Series:
    scores = {}

    sorted_meta = metadata.sort_values(by=["time"])
    sorted_ssgsea = ssgsea.loc[sorted_meta.index]
    TIME_LIMIT = np.quantile(sorted_meta["time"], rmst_time_limit)

    for split_by in ssgsea.columns:
        name, rmst_value, _ = fast_split_by_median(
            split_by,
            sorted_ssgsea[split_by],
            sorted_meta["time"].values,
            sorted_meta["event"].values,
            TIME_LIMIT,
        )

        scores[name] = rmst_value

    return pd.Series(scores, name=tcga)


"""
Code from Discovery-Science 2023: https://github.com/biolab/Discovery-Science-2023/blob/main/method/rmst_diff.py
"""


def fit_KM_sequence(sample_indicator: np.ndarray, event_indicator: np.ndarray):
    """
    Fits a sequence of survival probabilities where each step is unique timepoint from the data.
    The method assumes that samples are ordered by their survival time.

    Arguments
    ---------
    sample_indicator: np.ndarray (n_samples,)
        Indicator with 1 if sample is in cohort and 0 if it is not.
    event_indicator: np.ndarray (n_samples,)
        Indicator with 1 if event happened and 0 if it is censored.

    Returns
    -------
    np.ndarray (n_samples,)
        Survival probability of the cohort where indices correspond to the time points in the data.
    """
    n_series = np.cumsum(sample_indicator[::-1])[::-1]
    return np.append(
        [1], np.cumprod((n_series - event_indicator * sample_indicator) / n_series)
    )


def difference_RMST(km1: np.ndarray, km2: np.ndarray, time_values: np.ndarray):
    """
    Calculates the difference between the KM1 and KM1 restricted mean survival time.
    Basically integrates over the curves.
    KM curves are of equal length by design.

    Arguments
    ---------
    km1: np.ndarray (n_samples, )
        KM 1 curve with survival probabilities for each time point.
    km2: np.ndarray (n_samples, )
        KM 2 curve with survival probabilities for each time point.
    time_values: np.ndarray (n_samples, )
        Array of survival time points in the dataset

    Returns
    -------
    float
        The difference in RMST (AOC of difference in RMST)
    """
    LEN = len(time_values)  # the last value is the maximal (resticted at) time
    dt = np.ediff1d(time_values, to_begin=time_values[0])
    return np.sum((km1[:LEN] - km2[:LEN]) * dt)


def fast_split_by_median(
    feature, feature_values, sorted_time, sorted_events, TIME_LIMIT
):
    # Split values by median
    cutoff = feature_values.median()
    strata = (feature_values > cutoff).astype(bool)

    # Check if the highest time in strata is lower than time limit
    #   If so, RMST would be undefined or misleading at best.
    time_limit = min(sorted_time[strata].max(), sorted_time[~strata].max())

    limit_is_lower = True
    if TIME_LIMIT <= time_limit:
        time_limit = TIME_LIMIT
    else:
        print("TIME_LIMIT is lower than maximal time in one of the cohorts.")
        limit_is_lower = False

    # creates KM sequences for both strata
    km1 = fit_KM_sequence(strata, sorted_events)
    km2 = fit_KM_sequence(~strata, sorted_events)

    # calculates the difference between those
    dif_rstm = difference_RMST(
        km1, km2, list(sorted_time[sorted_time <= time_limit]) + [time_limit]
    )

    return (feature, dif_rstm, limit_is_lower)
