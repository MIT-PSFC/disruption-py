from disruption_py.database import *
from sklearn.model_selection import train_test_split

#TODO: Consider renaming 'target' to 'label' to be consistent with other code
DEFAULT_COLS = [
    'ip_error_frac',
    'Greenwald_fraction',
    'n_equal_1_normalized',
    'q95',
    'li',
    'radiated_fraction',
    'kappa',
    'dipprog_dt',
    'intentional_disruption',
    'power_supply_railed',
    'other_hardware_failure',
    'shot',
    'time_until_disrupt',
    'time']
REQUIRED_COLS = ['time', 'time_until_disrupt', 'shot']
TOKAMAKS = {'d3d': create_d3d_handler,
            'cmod': create_cmod_handler, 'east': create_east_handler}
DEFAULT_THRESHOLD =0.35 # Time until disrupt threshold for binary classification
BLACK_WINDOW_THRESHOLD = 5.e-3 # A 'black window' threshold [s]; obscures input data from a window in time on disruptive shots during trianing/testing
DEFAULT_RATIO = .2 # Ratio of test data to total data and validation data to train data

def create_target(time_until_disrupt, threshold, label_type, multiple_thresholds):
    if label_type == 'binary':
        dis = np.where((time_until_disrupt > threshold) & (~np.isnan(time_until_disrupt)))[0]  
        ndis = np.where(np.isnan(time_until_disrupt))[0] 
        target = np.ones(np.size(time_until_disrupt), dtype=int)
        target[np.hstack([dis, ndis])] = 0
    else:
        raise NotImplementedError('Only binary labels are implemented')


def get_dataset_df(cols=DEFAULT_COLS, shot_ids=None,threshold = DEFAULT_THRESHOLD tokamak='d3d', required_cols=REQUIRED_COLS):
    tokamak = TOKAMAKS[tokamak]()
    dataset_df = tokamak.get_shot_data(shot_ids, cols)
    cols = list(dataset_df.columns)
    if not set(required_cols).issubset(set(cols)):
        raise ValueError('Required columns not in dataset')
    # for col in required_cols:
    #     cols.remove(col)
    # cols = sorted(cols)
    # cols = cols + required_cols
    dataset_df['target'] = create_target(dataset_df['time_until_disrupt'].values, threshold, 'binary', False)
    

def filter_dataset_df(df, **kwargs):
    if 'exclude_non_disruptive' in kwargs and kwargs['exclude_non_disruptive']:
        keep_disrupt = np.where(~np.isnan(df['time_until_disrupt'].values))[0]
        df = df.iloc[keep_disrupt]
    if 'exclude_black_window' in kwargs and kwargs['exclude_black_window'] != 0:
        df = df[(df['time_until_disrupt'] > kwargs['exclude_black_window']) | np.isnan(df['time_until_disrupt'])]
    if 'impute' in kwargs and kwargs['impute']:
        raise NotImplementedError('Imputation not implemented')
    else:
        features = list(df.columns)
        df = df.dropna(subset=[feat for feat in features  if feat != 'time_until_disrupt'])
    if 'print_to_csv' in kwargs:
        if kwargs['print_to_csv']:
            df.to_csv('filtered_dataset.csv')
    return df

#TODO: Add other train_test_split options besides random 
def create_dataset(df, **kwargs):
    features = list(df.columns)
    features.remove('target')
    X, y = df.loc[features], df['target']
    ratio = kwargs['ratio'] if 'ratio' in kwargs else DEFAULT_RATIO
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=ratio, random_state=42)
    return X_train, X_test, y_train, y_test

def main():
    dataset_df = get_dataset_df()
    dataset_df = filter_dataset_df(dataset_df, exclude_non_disruptive=True, exclude_black_window=BLACK_WINDOW_THRESHOLD)
    X_train, X_test, y_train, y_test = create_dataset(dataset_df, ratio=DEFAULT_RATIO)
    print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
    df.to_csv('./whole_df.csv', sep=',', index=False)
    df_train = pd.concat([X_train, y_train], axis=1)
    df_train.to_csv('./train.csv', sep=',', index=False)
    df_test = pd.concat([X_test, y_test], axis=1)
    df_test.to_csv('./test.csv', sep=',', index=False)
    
