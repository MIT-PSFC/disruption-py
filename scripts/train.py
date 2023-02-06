import argparse

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix

def create_model(model_type):
    if model_type == 'random_forest':
        model = RandomForestClassifier(n_estimators=100, max_depth=2, random_state=0)
    else:
        raise NotImplementedError('Only random forest implemented')
    return model

def train_local(x_train, y_train,x_test,y_test, **kwargs):
    if 'features' in kwargs and kwargs['features'] is not None:
        features = kwargs['features']
        x_train = x_train[features]
    if 'model' in kwargs and kwargs['model'] is not None:
        raise NotImplementedError('Only random forest implemented')
    else:
        model = RandomForestClassifier(n_estimators=100, max_depth=2, random_state=0)
        model.fit(x_train, y_train)
        return {'model': model,
        'predictions': model.predict(x_train),
        'predictions_proba': model.predict_proba(x_train),
        'features': model.feature_importances_,
        'std_imp':np.std([model.feature_importances_ for model in forest.estimators_], axis=0),
        'score': model.score(x_test, y_test)}

def main(args):
    train = pd.read_csv(args.train_path)
    test = pd.read_csv(args.test_path)
    y_train = train['label']
    x_train = train.drop('label', axis=1)
    y_test = test['label']
    x_test = test.drop('label', axis=1)
    results_dict = train_local(x_train, y_train, x_test, y_test, features=args.features)
    conf_mat = confusion_matrix(y_test, results_dict['predictions'])
    print(conf_mat)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Train a model on a dataset.")
    parser.add_argument('--train_path', type=str, help='Path to training data. Must be a csv file.', default='./train.csv')
    parser.add_argument('--test_path', type=str, help='Path to test data. Must be a csv file.', default='./test.csv')
    parser.add_argument('--features', type=str, nargs='+', help='List of features to use for training. If not provided, all features will be used.', default=None)
    args = parser.parse_args()
    main(args)
