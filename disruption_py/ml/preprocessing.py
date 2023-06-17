import logging
import random
import os 

import numpy as np
import pandas as pd 
from sklearn.model_selection import train_test_split
from sklearn.impute import SimpleImputer

from disruption_py.database import create_d3d_handler, create_cmod_handler, create_east_handler
from disruption_py.shots import D3DShot, CModShot
from disruption_py.utils import generate_id, exp_filter

DEFAULT_COLS = [
    #    'ip_error_frac', Use just need to fix bugs
    'Greenwald_fraction',
    'n_equal_1_normalized',
    'q95',
    'li',
    'radiated_fraction',
    'kappa',
    # 'dipprog_dt',
    # 'intentional_disruption',
    # 'power_supply_railed',
    # 'other_hardware_failure',
    'shot',
    'time_until_disrupt',
    'time']
PAPER_COLS = [
    'aminor',
    'n_e',
    'ip',
    'delta',
    'li',
    'Wmhd',
    'kappa',
    # 'squareness'
]
DERIVED_PAPER_COLS = [
    'ip-exp-10-none',
    'ip-exp-50-none',
]
REQUIRED_COLS = ['time', 'time_until_disrupt', 'shot']
TOKAMAKS = {'d3d': create_d3d_handler,
            'cmod': create_cmod_handler, 'east': create_east_handler}
SHOT_CLASSES = {'d3d': D3DShot,
                'cmod': CModShot, 'east': None}
DEFAULT_THRESHOLD = 0.35  # Time until disrupt threshold for binary classification
# A 'black window' threshold [s]; obscures input data from a window in time on disruptive shots during trianing/testing
BLACK_WINDOW_THRESHOLD = 5.e-3
DEFAULT_RATIO = .2  # Ratio of test data to total data and validation data to train data

LOGGER = logging.getLogger('disruption_py')
def create_label(time_until_disrupt, threshold, label_type, multiple_thresholds):
    if label_type == 'binary':
        print(time_until_disrupt)
        dis = np.where((time_until_disrupt > threshold) &
                       (~np.isnan(time_until_disrupt)))[0]
        ndis = np.where(np.isnan(time_until_disrupt))[0]
        target = np.ones(np.size(time_until_disrupt), dtype=int)
        target[np.hstack([dis, ndis])] = 0
    else:
        raise NotImplementedError('Only binary labels are implemented')
    return target


# TODO: Add support for CMOD
def get_dataset_df(data_source=2, cols=DEFAULT_COLS, efit_tree=None, shot_ids=None, threshold=DEFAULT_THRESHOLD, tokamak='d3d', required_cols=REQUIRED_COLS, label='binary', **kwargs):
    if tokamak not in ['d3d', 'cmod']:
        raise NotImplementedError(
            "Currently only support DIII-D and Alcator C-MOD data retrieval")
    else:
        tokamak_handler = TOKAMAKS[tokamak]()
    timebase_signal = kwargs.get('timebase_signal', None)
    populate = kwargs.get('populate', 'default')
    label = kwargs.get('label', 'none')
    if shot_ids is None:
        random_state = kwargs.get('random_state', 8808)
        shot_ids = tokamak_handler.get_disruption_table_shotlist()['shot']
        random.Random(random_state).shuffle(shot_ids)
        shot_ids = shot_ids[:1200]
    if data_source == 0:
        shots = [tokamak_handler.get_shot(shot_id, efit_tree) for shot_id in shot_ids]
        dataset_df = pd.concat([shot.data for shot in shots])
    elif data_source == 1:
        raise NotImplementedError
    elif data_source == 2:
        dataset_df = tokamak_handler.get_shot_data(shot_ids, cols)
    elif data_source == 3:
        shots = []
        timebase_signal = kwargs.get('timebase_signal', None)
        for idx,shot_id in enumerate(shot_ids):
            percent_complete = idx/len(shot_ids)*100
            try:
                if tokamak == 'd3d':
                    if efit_tree is None:
                        shots.append(D3DShot(shot_id, tokamak_handler.get_efit_tree(
                            shot_id), disruption_time=tokamak_handler.get_disruption_time(shot_id), timebase_signal=timebase_signal, populate=populate))
                    else:
                        shots.append(D3DShot(shot_id, efit_tree, disruption_time=tokamak_handler.get_disruption_time(shot_id),
                                             timebase_signal=timebase_signal, populate=populate))
                elif tokamak == 'cmod':
                    #shots.append(CModShot("cmod",shot_id=shot_id))
                    shots.append(CModShot(shot_id=shot_id, disruption_time=tokamak_handler.get_disruption_time(shot_id)))
                LOGGER.info(f"[Shot {shot_id}]:Generated shot object, {idx} of {len(shot_ids)} ({percent_complete:.1f}% percent complete)' ")
            except Exception as e:
                LOGGER.info(f"[Shot {shot_id}]:Failed to generate shot object, {idx} of {len(shot_ids)} ({percent_complete:.1f}% percent complete)'")
                # exc_type, exc_obj, exc_tb = sys.exc_info()
                # fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                LOGGER.debug(f"[Shot {shot_id}]:{e}")
        dataset_df = pd.concat([shot.data for shot in shots])[cols]
    else:
        raise ValueError(
            'Datasource must be one of 4 options: 0,1,2,3. See generate_datsets.py -h for more details.')
    LOGGER.info(f"Successfully processed {len(shots)}/{len(shot_ids)} shots")
    dataset_df = dataset_df.fillna(value=np.nan)
    cols = list(dataset_df.columns)
    if not set(required_cols).issubset(set(cols)):
        raise ValueError('Required columns not in dataset')
    if label != 'none':
        dataset_df['label'] = create_label(
            dataset_df['time_until_disrupt'].values, threshold, label, False)
    else:
        dataset_df['label'] = np.nan
    return dataset_df


def add_derived_features(df, cols):
    for col in cols:
        feature, func, *args = col.split('-')
        if func == 'exp':
            w = float(args[0])/100.
            df[col] = exp_filter(df[feature], w, *args[1:])
        else:
            print(f"{func} is not supported")
    return df


def filter_dataset_df(df, **kwargs):
    exclude_non_disrupted = kwargs.get('exclude_non_disruptive', False)
    exclude_black_window = kwargs.get('exclude_black_window', 0)
    impute = kwargs.get('impute', True)
    write_to_csv = kwargs.get('write_to_csv', False)
    csv_path = kwargs.get('csv_path', r'./output/filtered_dataset.csv')
    if exclude_non_disrupted:
        keep_disrupt = np.where(~np.isnan(df['time_until_disrupt'].values))[0]
        df = df.iloc[keep_disrupt]
    if exclude_black_window != 0:
        df = df[(df['time_until_disrupt'] > exclude_black_window)
                | np.isnan(df['time_until_disrupt'])]
    print(df.head())
    features = list(df.columns)
    if impute:
        df = impute_shot_df_NaNs(df)
        df = df.dropna(
            subset=[feat for feat in features if feat != 'time_until_disrupt'])
    else:
        df = df.dropna(
            subset=[feat for feat in features if feat != 'time_until_disrupt'])
    if write_to_csv:
        df.to_csv(csv_path)
    return df


def create_dataset(df, split_by_shot=True, **kwargs):
    features = list(df.columns)
    features.remove('label')
    X, y = df[features], df['label']
    ratio = kwargs.get('ratio', DEFAULT_RATIO)
    random_state = kwargs.get('random_state', 42)
    if split_by_shot:
        shots = pd.unique(df['shot'])
        random.Random(random_state).shuffle(shots)
        n_test = int(np.ceil(len(shots)*ratio))
        test_shots = shots[:n_test]
        train_shots = shots[n_test:]
        X_train, X_test = X[df['shot'].isin(train_shots)], X[df['shot'].isin(
            test_shots)]
        y_train, y_test = y[df['shot'].isin(train_shots)], y[df['shot'].isin(
            test_shots)]
    else:
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=ratio, random_state=random_state)
    return X_train, X_test, y_train, y_test


def parse_feature_cols(feature_str):
    if feature_str is None:
        return PAPER_COLS, []
    elif os.path.isfile(feature_str):
        all_cols = pd.read_csv(
            feature_str, header=None).iloc[:, 0].values
    else:
        all_cols = [feature.strip()
                    for feature in feature_str.split(",")]
    feature_cols = []
    derived_feature_cols = []
    for col in all_cols:
        if '-' in col:
            derived_feature_cols.append(col)
        else:
            feature_cols.append(col)
    return feature_cols, derived_feature_cols

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