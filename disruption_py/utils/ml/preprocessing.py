import os
import random
import pandas as pd
import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.model_selection import train_test_split
from disruption_py.utils.constants import DEFAULT_RATIO, DEFAULT_THRESHOLD, PAPER_COLS
from disruption_py.utils.math_utils import exp_filter


def parse_feature_cols(feature_str):
    if feature_str is None:
        return PAPER_COLS, []
    elif os.path.isfile(feature_str):
        all_cols = pd.read_csv(feature_str, header=None).iloc[:, 0].values
    else:
        all_cols = [feature.strip() for feature in feature_str.split(",")]
    feature_cols = []
    derived_feature_cols = []
    for col in all_cols:
        if "-" in col:
            derived_feature_cols.append(col)
        else:
            feature_cols.append(col)
    return feature_cols, derived_feature_cols


def create_label(
    time_until_disrupt,
    threshold=DEFAULT_THRESHOLD,
    label_type="binary",
    multiple_thresholds=None,
):
    if label_type == "binary":
        dis = np.where((time_until_disrupt > threshold) | np.isnan(time_until_disrupt))[
            0
        ]
        # ndis = np.where(np.isnan(time_until_disrupt))[0]
        target = np.ones(np.size(time_until_disrupt), dtype=int)
        target[dis] = 0
    else:
        raise NotImplementedError("Only binary labels are implemented")
    return target


def add_derived_features(df, cols):
    for col in cols:
        feature, func, *args = col.split("-")
        if func == "exp":
            w = float(args[0]) / 100.0
            df = df.sort_values(["shot", "time"])
            df[col] = df.groupby("shot")[feature].transform(
                lambda x: exp_filter(x.reset_index(drop=True), w, *args[1:])
            )
        else:
            print(f"{func} is not supported")
    return df


def create_dataset(df, split_by_shot=True, **kwargs):
    features = list(df.columns)
    features.remove("label")
    X, y = df[features], df["label"]
    ratio = kwargs.get("ratio", DEFAULT_RATIO)
    random_state = kwargs.get("random_state", 42)
    if split_by_shot:
        shots = pd.unique(df["shot"])
        random.Random(random_state).shuffle(shots)
        n_test = int(np.ceil(len(shots) * ratio))
        test_shots = shots[:n_test]
        train_shots = shots[n_test:]
        X_train, X_test = (
            X[df["shot"].isin(train_shots)],
            X[df["shot"].isin(test_shots)],
        )
        y_train, y_test = (
            y[df["shot"].isin(train_shots)],
            y[df["shot"].isin(test_shots)],
        )
    else:
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=ratio, random_state=random_state
        )
    return X_train, X_test, y_train, y_test


def filter_dataset_df(df, **kwargs):
    exclude_non_disrupted = kwargs.get("exclude_non_disruptive", False)
    exclude_black_window = kwargs.get("exclude_black_window", 0)
    impute = kwargs.get("impute", True)
    write_to_csv = kwargs.get("write_to_csv", False)
    csv_path = kwargs.get("csv_path", r"./output/filtered_dataset.csv")
    if exclude_non_disrupted:
        keep_disrupt = np.where(~np.isnan(df["time_until_disrupt"].values))[0]
        df = df.iloc[keep_disrupt]
    if exclude_black_window != 0:
        df = df[
            (df["time_until_disrupt"] > exclude_black_window)
            | np.isnan(df["time_until_disrupt"])
        ]
    features = list(df.columns)
    if impute:
        df = impute_shot_df_NaNs(df)
        df = df.dropna(
            subset=[feat for feat in features if feat != "time_until_disrupt"]
        )
    else:
        df = df.dropna(
            subset=[feat for feat in features if feat != "time_until_disrupt"]
        )
    if write_to_csv:
        df.to_csv(csv_path)
    return df


# TODO: Clean up this function
# TODO: Fix issue with None
def impute_shot_df_NaNs(df, strategy="median", missing_values=np.nan, cutoff=0.5):
    """
    This routine imputes missing values for all columns in a given df.
    Values are imputed on the basis of different shot numbers.
    It returns a dataframe where all columns that have missing values
    now have them imputed according to the selected strategy. In addition,
    for every imputed column, a new flag_ column is added: if an original
    NaN was in that row, a 1 is there otherwise 0.
    df:                pandas dataframe;
    strategy:          kwarg, default is 'median'. Other possible ones are 'mean' or 'nearest';
    missing_values:    kwarg, default is 'NaN';
    return df
    """

    impute = SimpleImputer(missing_values=missing_values, strategy=strategy)
    variables = []
    ordered_names = list(df.columns)
    cols = [col for col in ordered_names if col not in ["time_until_disrupt", "label"]]
    for col in cols:
        if df[col].isnull().values.any():
            variables.append(col)

    shots = np.unique(df["shot"].values)

    for var in variables:
        imputed_var = []
        flag_imputed_nans = []

        for shot in shots:
            index = np.where(df["shot"].values == shot)[0]
            original_var = np.array(df[var].values[index], dtype=np.float64)
            tmp_flag = np.zeros(np.size(original_var))
            nan_original_var = np.where(np.isnan(original_var))[0]
            # if the whole column is NaN or more than 50% are NaNs
            # then the whole shot should be discarded
            if (np.size(nan_original_var) == np.size(index)) or (
                np.size(nan_original_var) >= cutoff * np.size(index)
            ):
                tmp_var = original_var * np.nan
                tmp_flag = original_var * np.nan
            else:
                tmp_var = impute.fit_transform(original_var.reshape(-1, 1))
                tmp_flag[nan_original_var] = 1

            imputed_var.extend(np.hstack(tmp_var))
            flag_imputed_nans.extend(tmp_flag)

        df.loc[:, var] = imputed_var
    return df
