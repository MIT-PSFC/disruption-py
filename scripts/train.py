import argparse
import pickle

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, fbeta_score

def create_model(model_type):
    if model_type == 'random_forest':
        model = RandomForestClassifier(n_estimators=100, max_depth=10, random_state=0)
    else:
        raise NotImplementedError('Only random forest implemented')
    return model

def train_local(x_train, y_train,x_test,y_test, **kwargs):
    if 'omit_after_quench' in kwargs and kwargs['omit_after_quench'] is not None:
        raise NotImplementedError("Talk to Cristina about implementing")
    if 'features' in kwargs and kwargs['features'] is not None:
        features = kwargs['features']
        x_train = x_train[features]
    if 'model' in kwargs and kwargs['model'] is not None:
        raise NotImplementedError('Only random forest implemented')
    else:
        model = RandomForestClassifier(n_estimators=500, max_depth=10, random_state=0)

        model.fit(x_train, y_train)
        return {'model': model,
        'predictions': model.predict(x_train),
        'predictions_proba': model.predict_proba(x_train),
        'features': model.feature_importances_,
        'std_imp':np.std([tree.feature_importances_ for tree in model.estimators_], axis=0),
        'score': model.score(x_test, y_test)}

def main(args):
    if 'pkl' in args.train_path:
        with open(args.train_path, 'rb') as f:
            train = pickle.load(f,encoding='latin1')
        x_train = train['X_train']
        y_train = train['y_train'].astype(int)
    else:
        train = pd.read_csv(args.train_path)
        y_train = train['label']
        x_train = train.drop('label', axis=1)
    if 'pkl' in args.test_path:
        with open(args.test_path,'rb') as f:
            test = pickle.load(f,encoding='latin1')
        x_test = test['X_test']
        y_test = test['y_test']
    else:
        test = pd.read_csv(args.test_path)
        y_test = test['label']
        x_test = test.drop('label', axis=1)
    results_dict = train_local(x_train, y_train, x_test, y_test, features=args.features)
    model = results_dict['model']  
    test_predictions = model.predict(x_test)  
    print("Train F2-Score:", fbeta_score(y_train, results_dict['predictions'], beta=2))
    print("Test F2-Score:", fbeta_score(y_test, test_predictions, beta=2))
    conf_mat = confusion_matrix(y_train, results_dict['predictions'])
    disp_train = ConfusionMatrixDisplay(confusion_matrix = conf_mat, display_labels = model.classes_)
    conf_mat = confusion_matrix(y_test,test_predictions)
    disp_test = ConfusionMatrixDisplay(confusion_matrix=conf_mat, display_labels=model.classes_)
    disp_train.plot()
    disp_test.plot()
    disp_train.ax_.set_title('Train Matrix')
    disp_test.ax_.set_title('Test Matrix')
    plt.show()
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Train a model on a dataset.")
    parser.add_argument('--train_path', type=str, help='Path to training data. Must be a csv file.', default='./train.csv')
    parser.add_argument('--test_path', type=str, help='Path to test data. Must be a csv file.', default='./test.csv')
    parser.add_argument('--features', type=str, nargs='+', help='List of features to use for training. If not provided, all features will be used.', default=None)
    args = parser.parse_args()
    main(args)
