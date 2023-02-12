import random
import argparse
try:
    import importlib.resources as importlib_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources

from sklearn.model_selection import train_test_split

from disruption_py.database import *
from disruption_py.utils import impute_shot_df_NaNs, exp_filter
import disruption_py.data

# TODO: Add argument for picking between MDSPlus data and SQL data
# TODO: Add argument for signal to use as master timebase
# TODO  Add argument for efit to use (get_efit_tree(shot_id) or a specific efit tree name)
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
    'squareness'
]
DERIVED_PAPER_COLS = [
    'ip-exp-10-none',
    'ip-exp-50-none',
]
REQUIRED_COLS = ['time', 'time_until_disrupt', 'shot']
TOKAMAKS = {'d3d': create_d3d_handler,
            'cmod': create_cmod_handler, 'east': create_east_handler}
DEFAULT_THRESHOLD = 0.35  # Time until disrupt threshold for binary classification
# A 'black window' threshold [s]; obscures input data from a window in time on disruptive shots during trianing/testing
BLACK_WINDOW_THRESHOLD = 5.e-3
DEFAULT_RATIO = .2  # Ratio of test data to total data and validation data to train data


def create_label(time_until_disrupt, threshold, label_type, multiple_thresholds):
    if label_type == 'binary':
        dis = np.where((time_until_disrupt > threshold) &
                       (~np.isnan(time_until_disrupt)))[0]
        ndis = np.where(np.isnan(time_until_disrupt))[0]
        target = np.ones(np.size(time_until_disrupt), dtype=int)
        target[np.hstack([dis, ndis])] = 0
    else:
        raise NotImplementedError('Only binary labels are implemented')
    return target

# TODO: Add support for CMOD


def get_dataset_df(data_source=2, cols=DEFAULT_COLS, efit_tree=None, shot_ids=None, threshold=DEFAULT_THRESHOLD, tokamak='d3d', required_cols=REQUIRED_COLS):
    if tokamak != 'd3d':
        raise NotImplementedError(
            "Currently only support DIII-D data retrieval")
    tokamak = TOKAMAKS[tokamak]()
    if shot_ids is None:
        shot_ids = tokamak.get_disruption_table_shotlist()['shot']
        random.shuffle(shot_ids)
        shot_ids = shot_ids[:1200]
    if data_source == 0:
        shots = [tokamak.get_shot(shot_id, efit_tree) for shot_id in shot_ids]
        dataset_df = pd.concat([shot.data for shot in shots])
    elif data_source == 1:
        raise NotImplementedError
    elif data_source == 2:
        dataset_df = tokamak.get_shot_data(shot_ids, cols)
    elif data_source == 3:
        shots = []
        for shot_id in shot_ids:
            if efit_tree is None:
                shots.append(D3DShot(shot_id, tokamak.get_efit_tree(shot_id)))
            else:
                shots.append(D3DShot(shot_id, efit_tree))
        dataset_df = pd.concat([shot.data for shot in shots])
    else:
        raise ValueError(
            'Datasource must be one of 4 options: 0,1,2,3. See generate_datsets.py -h for more details.')
    dataset_df = dataset_df.fillna(value=np.nan)
    cols = list(dataset_df.columns)
    if not set(required_cols).issubset(set(cols)):
        raise ValueError('Required columns not in dataset')
    # for col in required_cols:
    #     cols.remove(col)
    # cols = sorted(cols)
    # cols = cols + required_cols
    dataset_df['label'] = create_label(
        dataset_df['time_until_disrupt'].values, threshold, 'binary', False)
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
    if 'exclude_non_disruptive' in kwargs and kwargs['exclude_non_disruptive']:
        keep_disrupt = np.where(~np.isnan(df['time_until_disrupt'].values))[0]
        df = df.iloc[keep_disrupt]
    if 'exclude_black_window' in kwargs and kwargs['exclude_black_window'] != 0:
        df = df[(df['time_until_disrupt'] > kwargs['exclude_black_window'])
                | np.isnan(df['time_until_disrupt'])]
    if 'impute' in kwargs and kwargs['impute']:
        df = impute_shot_df_NaNs(df)
        df = df.dropna()
    else:
        features = list(df.columns)
        df = df.dropna(
            subset=[feat for feat in features if feat != 'time_until_disrupt'])
    if 'print_to_csv' in kwargs:
        if kwargs['print_to_csv']:
            df.to_csv('filtered_dataset.csv')
    return df


def create_dataset(df, split_by_shot=True, **kwargs):
    features = list(df.columns)
    features.remove('label')
    features.remove('time_until_disrupt')
    features.remove('shot')
    X, y = df[features], df['label']
    ratio = kwargs.get('ratio', DEFAULT_RATIO)
    if split_by_shot:
        shots = pd.unique(df['shot'])
        random.shuffle(shots)
        n_test = int(np.ceil(len(shots)*ratio))
        test_shots = shots[:n_test]
        train_shots = shots[n_test:]
        X_train, X_test = X[df['shot'].isin(train_shots)], X[df['shot'].isin(
            test_shots)]
        y_train, y_test = y[df['shot'].isin(train_shots)], y[df['shot'].isin(
            test_shots)]
    else:
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=ratio, random_state=42)
    return X_train, X_test, y_train, y_test


def parse_feature_cols(feature_str):
    if feature_str is None:
        return PAPER_COLS, DERIVED_PAPER_COLS
    elif os.is_file(feature_str):
        all_cols = pd.read_csv(
            feature_str, header=None).iloc[:, 0].values
    else:
        all_cols = [feature.trim()
                    for feature in feature_str.split(",")]
    feature_cols = []
    derived_feature_cols = []
    for col in all_cols:
        if '-' in col:
            derived_feature_cols.append(col)
        else:
            feature_cols.append(col)
    return feature_cols, derived_feature_cols


def main(args):
    if args.shotlist is None:
        with importlib_resources.path(disruption_py.data, "paper_shotlist.txt") as p:
            args.shotlist = str(p)
    shot_ids = pd.read_csv(args.shotlist, header=None).iloc[:, 0].values
    feature_cols, derived_feature_cols = parse_feature_cols(args.feature_cols)
    dataset_df = get_dataset_df(
        feature_cols + REQUIRED_COLS, shot_ids=shot_ids)
    dataset_df = add_derived_features(dataset_df, derived_feature_cols)
    dataset_df = filter_dataset_df(dataset_df, exclude_non_disruptive=True,
                                   exclude_black_window=BLACK_WINDOW_THRESHOLD, impute=True)
    X_train, X_test, y_train, y_test = create_dataset(
        dataset_df, ratio=DEFAULT_RATIO)
    print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
    dataset_df.to_csv('./whole_df.csv', sep=',', index=False)
    df_train = pd.concat([X_train, y_train], axis=1)
    df_train.to_csv('./train.csv', sep=',', index=False)
    df_test = pd.concat([X_test, y_test], axis=1)
    df_test.to_csv('./test.csv', sep=',', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate DPRF compatible datasets for training and inference. Currently only supports DIII-D data")
    parser.add_argument('--shotlist', type='str',
                        help='Path to file specifying shotlist', default=None)
    parser.add_argument('--feature_cols', type='str',
                        help='Either a file or comma-separated list of desired feature columns', default=None)
    parser.add_argument('--output_dir', type=str,
                        help='Path to generated data.', default='./data/')
    parser.add_argument('--timebase_signal', type=str,
                        help='Signal whose timebase will be used as the unifying timebase of the dataset.', default=None)
    parser.add_argument('--efit_tree', type=str,
                        help="Name of efit tree to use for each shot. If left as None, the script will use the get_efit_tree method in database.py.", default=None)
    parser.add_argument('--data_source', type=int, choices=[
                        0, 1, 2, 3], help=r"0: Default to SQL database then MDSPlus.\n1: Default to MDSPlus then SQL database.\n2: SQL database only.\n3: MDSPlus only.", default=0)
    args = parser.parse_args()
    main(args)
