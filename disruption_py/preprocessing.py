import numpy as np
from sklearn.impute import SimpleImputer

# TODO: Clean up this function
# TODO: Fix issue with None


def impute_shot_df_NaNs(df, strategy='median', missing_values=np.nan):
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

    impute = SimpleImputer(missing_values=missing_values,
                           strategy=strategy, verbose=0)
    variables = []
    ordered_names = list(df.columns)
    cols = [col for col in ordered_names if col not in [
        'time_until_disrupt', 'label']]
    for col in cols:
        if df[col].isnull().values.any():
            variables.append(col)

    shots = np.unique(df['shot'].values)

    for var in variables:
        imputed_var = []
        flag_imputed_nans = []

        for shot in shots:
            index = np.where(df['shot'].values == shot)[0]
            original_var = np.array(df[var].values[index], dtype=np.float64)
            tmp_flag = np.zeros(np.size(original_var))
            nan_original_var = np.where(np.isnan(original_var))[0]
            # if the whole column is NaN or more than 50% are NaNs
            # then the whole shot should be discarded
            if (np.size(nan_original_var) == np.size(index)) or (np.size(nan_original_var) >= 0.8 * np.size(index)):
                tmp_var = original_var * np.nan
                tmp_flag = original_var * np.nan
            else:
                tmp_var = impute.fit_transform(original_var.reshape(-1, 1))
                tmp_flag[nan_original_var] = 1

            imputed_var.extend(np.hstack(tmp_var))
            flag_imputed_nans.extend(tmp_flag)

        df.loc[:, var] = imputed_var
    return df
