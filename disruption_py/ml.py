from disruption_py.database import *
from disruption_py.utils import impute_shot_df_NaNs, exp_filter, generate_id

import random
import os
import pickle
import joblib
try:
    import importlib.resources as importlib_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources
from datetime import date

import pandas as pd
import numpy as np
import h5py
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier


DEFAULT_ORDER = {
    'ip': 0,
    'Wmhd': 1,
    'n_e': 2,
    'aminor': 3,
    'delta': 4,
    'squareness': 5,
    'kappa': 6,
    'li': 7,
}

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

LOGGER = logging.getLogger('disruption_py')

###############
# Generating #
#  Datasets  #
###############


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
    if tokamak != 'd3d':
        raise NotImplementedError(
            "Currently only support DIII-D data retrieval")
    tokamak = TOKAMAKS[tokamak]()
    timebase_signal = kwargs.get('timebase_signal', None)
    populate = kwargs.get('populate', 'default')
    label = kwargs.get('label', 'none')
    if shot_ids is None:
        random_state = kwargs.get('random_state', 8808)
        shot_ids = tokamak.get_disruption_table_shotlist()['shot']
        random.Random(random_state).shuffle(shot_ids)
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
        timebase_signal = kwargs.get('timebase_signal', None)
        for shot_id in shot_ids:
            try:
                if efit_tree is None:
                    shots.append(D3DShot(shot_id, tokamak.get_efit_tree(
                        shot_id), disruption_time=tokamak.get_disruption_time(shot_id), timebase_signal=timebase_signal, populate=populate))
                else:
                    shots.append(D3DShot(shot_id, efit_tree, disruption_time=tokamak.get_disruption_time(shot_id),
                                         timebase_signal=timebase_signal, populate=populate))
                LOGGER.info(f"[Shot {shot_id}]:Generated shot object")
            except Exception as e:
                LOGGER.info(f"[Shot {shot_id}]:Failed to generate shot object")
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


def get_d3d_dataset_df(data_source=2, cols=DEFAULT_COLS, efit_tree=None, shot_ids=None, threshold=DEFAULT_THRESHOLD, tokamak='d3d', required_cols=REQUIRED_COLS, label='binary', **kwargs):


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

###############
#   Models    #
###############


class HDF5RFModel():
    def __init__(self, n_trees, n_nodes, n_classes, value_arr, tree_start_arr, feature_arr, children_left_arr, children_right_arr, threshold_arr, stride):
        self.n_trees = n_trees
        self.n_nodes = n_nodes
        self.n_classes = n_classes
        self.value_arr = value_arr
        self.tree_start_arr = tree_start_arr
        self.feature_arr = feature_arr
        self.children_left_arr = children_left_arr
        self.children_right_arr = children_right_arr
        self.threshold_arr = threshold_arr
        self.stride = stride
        self.classes_ = [0, 1]

    def predict_proba(self, X):
        predictions = np.zeros((len(X)//self.stride,))
        results = np.zeros((self.n_trees,))
        for j in np.arange(len(X)//self.stride):
            X_eval = X[self.stride*j, :]
            for i in np.arange(self.n_trees):
                current_node = 0
                current_index = self.tree_start_arr[i]+current_node
                while self.children_left_arr[current_index] != self.children_right_arr[current_index]:
                    if X_eval[self.feature_arr[current_index]] <= self.threshold_arr[current_index]:
                        current_node = self.children_left_arr[current_index]
                    else:
                        current_node = self.children_right_arr[current_index]
                    current_index = self.tree_start_arr[i]+current_node
                results[i] = self.value_arr[current_index][1] / \
                    np.sum(self.value_arr[current_index])
            predictions[j] = np.sum(results)/self.n_trees
        return predictions


def load_model_from_hdf5(model_path, model_type='RandomForestClassifier'):
    if model_type == 'RandomForestClassifier':
        f_dict = h5py.File(model_path, 'r')
        n_trees = f_dict['n_trees'][0]
        n_nodes = f_dict['n_nodes'][0]
        n_classes = f_dict['n_classes'][0]
        value_arr = f_dict['value'][:]
        tree_start_arr = f_dict['tree_start'][:]
        feature_arr = f_dict['feature'][:]
        children_left_arr = f_dict['children_left'][:]
        children_right_arr = f_dict['children_right'][:]
        threshold_arr = f_dict['threshold'][:]
        stride = 1
        return HDF5RFModel(n_trees, n_nodes, n_classes, value_arr, tree_start_arr,
                           feature_arr, children_left_arr, children_right_arr, threshold_arr, stride)
    else:
        raise ValueError('Model type not supported')


def create_model(model_type):
    if model_type == 'random_forest':
        model = RandomForestClassifier(
            n_estimators=245, max_depth=15, random_state=0)
    else:
        raise NotImplementedError('Only random forest implemented')
    return model


def load_model(model_path):
    if model_path.endswith('.joblib'):
        model = joblib.load(model_path)
    elif model_path.endswith('.pkl'):
        with open(model_path, 'rb') as f:
            model = pickle.load(f)
    elif model_path.endswith('.h5'):
        model = load_model_from_hdf5(model_path)
    else:
        raise ValueError(f"Unknown model file type: {model_path}")
    return model


def predict(model, df, order=DEFAULT_ORDER):
    arr = np.empty((len(df), len(order)))
    for col in order:
        arr[:, order[col]] = df[col].values
    predictions = model.predict_proba(arr)
    if len(predictions.shape) > 1:
        predictions = predictions[:, 1]
    df['score'] = predictions
    return df

#############
#   Train   #
#############


def grid_search(x_train, y_train, x_test, y_test, model_type, grid, **kwargs):
    """ Run a grid search over the given parameters.

    Parameters
    ----
    x_train: np.ndarray
        Training data
    y_train: np.ndarray
        Training labels
    x_test: np.ndarray
        Test data
    y_test: np.ndarray
        Test labels
    model_type: str
        Type of model to use
    grid: dict
        Dictionary of parameters to search over

    Returns
    ----
    list
        List of dictionaries with the parameters and results for each search
    """
    # Grid is a dictionary where each key is a parameter and each value is a list of values to try
    # kwargs are any other parameters to pass to the model
    searches = []
    for param in grid:
        for val in grid[param]:
            if len(searches) == 0:
                searches.append({param: val})
            else:
                searches = [dict(search, **{param: val})
                            for search in searches]
    results = []
    for search in searches:
        result = train_local(x_train, y_train, x_test,
                             y_test,  **search, **kwargs)
        results.append({'params': search, 'result': result})
    return results


def train_local(x_train, y_train, x_test, y_test, **kwargs):
    if 'omit_after_quench' in kwargs and kwargs['omit_after_quench'] is not None:
        raise NotImplementedError("Talk to Cristina about implementing")
    if 'features' in kwargs and kwargs['features'] is not None:
        features = kwargs['features']
        x_train = x_train[features]
    if 'model' in kwargs and kwargs['model'] is not None:
        raise NotImplementedError('Only random forest implemented')
    else:
        n_estimators = kwargs.pop('n_estimators', 245)
        max_depth = kwargs.pop('max_depth', 15)
        random_state = kwargs.pop('random_state', 0)
        n_jobs = kwargs.pop('n_jobs', -1)
        model = RandomForestClassifier(
            n_estimators=n_estimators, max_depth=max_depth, random_state=random_state, n_jobs=n_jobs, **kwargs)
        model.fit(x_train, y_train)
        return {'model': model,
                'predictions': model.predict(x_train),
                'predictions_proba': model.predict_proba(x_train),
                'features': model.feature_importances_,
                'std_imp': np.std([tree.feature_importances_ for tree in model.estimators_], axis=0),
                'score': model.score(x_test, y_test)}

#############
#   Eval    #
#############


def eval_shots(df, lower_threshold=.05, disruptivity=.45, window=.025):
    good_warnings = []
    missed_warnings = []
    false_alarms = []
    warning_times = dict()
    disruptions = []
    non_disruptions = []
    for _, shot_data in df.groupby(by=['shot']):
        shot_id = shot_data['shot'].values[0]
        time_until_disrupt = shot_data['time_until_disrupt'].values
        alarm = trigger_alarm(
            shot_data['time'].values, shot_data['score'].values, lower_threshold, disruptivity, window)
        plt.figure()
        plt.title(
            f"Shot {shot_id}\nlower_threshold={lower_threshold};disruptivity={disruptivity};window={window}")
        plt.plot(shot_data['time'], shot_data['score'],
                 label='Disruptivity Score')
        plt.plot(shot_data['time'], alarm, label='Alarm Activated')
        trigger = np.where(np.array(alarm) == 1)[0]
        disrupted = np.isnan(time_until_disrupt).all()
        if disrupted:
            disruptions.append(shot_id)
            plt.axvline(x=shot_data['time'].values[0] +
                        shot_data['time_until_disrupt'].values[0], label='Time of disruption')
        else:
            non_disruptions.append(shot_id)
        plt.legend()
        if np.size(trigger) > 0:
            if disrupted:
                warning_time = time_until_disrupt[trigger][0]
                warning_times[shot_id] = warning_time
                good_warnings.append(warning_time)
            else:
                false_alarms.append(shot_id)
        else:
            if disrupted:
                missed_warnings.append(shot_id)
    results = {}
    # TODO: Replace with f"" strings
    if len(disruptions) != 0:
        print(
            "Disruptions warned: %d/%d (%d%%)"
            % (len(good_warnings), len(disruptions), round(float(len(good_warnings)) / len(disruptions) * 100, 2))
        )
        print(
            "Missed disruptions: %d/%d (%d%%)"
            % (len(missed_warnings), len(disruptions), round(float(len(missed_warnings)) / len(disruptions) * 100, 2))
        )
    if len(non_disruptions) != 0:
        print(
            "False Alarms: %d/%d (%d%%)"
            % (len(false_alarms), len(non_disruptions), round(float(len(false_alarms)) / len(non_disruptions) * 100, 2))
        )
    results['good_warnings'] = good_warnings
    results['missed_warnings'] = missed_warnings
    results['false_alarms'] = false_alarms
    results['warning_times'] = warning_times
    results['TP'] = 0
    results['FN'] = 0
    results['FP'] = 0
    results['TN'] = 0
    if len(disruptions) > 0:
        results['TP'] = len(good_warnings) / float(len(disruptions))
        results['FN'] = len(missed_warnings) / float(len(disruptions))
    if len(non_disruptions) > 0:
        results['FP'] = len(false_alarms) / float(len(non_disruptions))
        results['TN'] = -(len(good_warnings) + len(false_alarms) +
                          len(missed_warnings) - len(df) / float(len(non_disruptions)))
    return results


def hist_time_th(thL, thH, thA, dt, v, v_old, c_old, a_old):
    # =================
    # double threshold
    # =================
    # if signal is below LOW threshold output is always 0
    if v < thL:
        v = 0
    # if signal is above HIGHT threshold output is always 1
    elif v > thH:
        v = 1
    else:
        v = v_old
    # =================
    # time accumulation
    # =================
    # when output signal is 0 counter gets reset
    if v == 0:
        c = 0
    # when output signal is 1 counter gets incremented
    else:
        c = c_old + dt
    # =================
    # sound alarm
    # =================
    # if counter surpasses allarm threshold or old allarm was on
    if c > thA or a_old:
        a = 1
    else:
        a = 0
    return v, c, a


def trigger_alarm(time, value, low_thr, high_thr, thA):
    output = []
    alarm = []
    if time[0] < 0:
        tstart = 0
    else:
        tstart = time[0]  # 0
    v_old = 0
    c_old = 0
    a_old = 0

    for t, v in zip(time, value):
        dt = t - tstart
        v, c, a = hist_time_th(low_thr, high_thr, thA,
                               dt, v, v_old, c_old, a_old)
        output.append(v)
        alarm.append(a)
        if output[0] != 0:
            v, c, a = 0, 0, 0
            output[0], alarm[0] = 0, 0
        v_old, c_old, a_old = v, c, a
        tstart = t
    return alarm
