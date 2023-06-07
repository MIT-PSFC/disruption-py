import argparse
import pickle
import joblib
import json
from datetime import date

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, fbeta_score

from disruption_py.utils import generate_id
from disruption_py.ml.train import grid_search, train_local

def main(args):
    if 'pkl' in args.train_path:
        with open(args.train_path, 'rb') as f:
            train = pickle.load(f, encoding='latin1')
        x_train = train['X_train']
        y_train = train['y_train'].astype(int)
    else:
        train = pd.read_csv(args.train_path)
        y_train = train['label']
        x_train = train.drop(
            ['label', 'time_until_disrupt', 'shot'], axis=1, errors='ignore')
    if 'pkl' in args.test_path:
        with open(args.test_path, 'rb') as f:
            test = pickle.load(f, encoding='latin1')
        x_test = test['X_test']
        y_test = test['y_test']
    else:
        test = pd.read_csv(args.test_path)
        y_test = test['label']
        x_test = test.drop(['label', 'time_until_disrupt',
                           'shot'], axis=1, errors='ignore')
    if args.grid_search is not None:
        grid = json.loads(args.grid_search)
        results = grid_search(x_train, y_train, x_test, y_test, **grid)
        with open(args.output_dir + f"grid_search_{args.unique_id}_.json", "w") as f:
            json.dump(results, f)
        return
    else:
        results = [train_local(
            x_train, y_train, x_test, y_test, features=args.features)]
        best_result = results[0]
    model = best_result['model']
    train['scores'] = best_result['predictions_proba']
    test['scores'] = model.predict(x_test)
    print("Train F2-Score:", fbeta_score(y_train,
          best_result['predictions'], beta=2))
    print("Test F2-Score:", fbeta_score(y_test, test['scores'], beta=2))
    conf_mat = confusion_matrix(y_train, best_result['predictions'])
    disp_train = ConfusionMatrixDisplay(
        confusion_matrix=conf_mat, display_labels=model.classes_)
    conf_mat = confusion_matrix(y_test, test['scores'])
    disp_test = ConfusionMatrixDisplay(
        confusion_matrix=conf_mat, display_labels=model.classes_)
    disp_train.plot()
    disp_test.plot()
    disp_train.ax_.set_title('Train Matrix')
    disp_test.ax_.set_title('Test Matrix')
    plt.show()
    joblib.dump(model, args.output_dir +
                f"random_forest{date.today()}_{args.unique_id}.joblib")
    with open(args.output_dir + f"train_{args.unique_id}_.json", "w") as f:
        json.dump(vars(args), f)
    print(f"Unique ID for this run: {args.unique_id}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Train a model on a dataset.")
    parser.add_argument('--train_path', type=str,
                        help='Path to training data. Must be a csv file.', default='./output/train.csv')
    parser.add_argument('--test_path', type=str,
                        help='Path to test data. Must be a csv file.', default='./output/test.csv')
    parser.add_argument('--output_dir', type=str, default='./output/',
                        help='Path to output model files. Must be a directory.')
    parser.add_argument('--grid_search', type=str, default=None,
                        help='Path to grid search parameters. Must be a json file.')
    parser.add_argument('--features', type=str, nargs='+',
                        help='List of features to use for training. If not provided, all features will be used.', default=None)
    parser.add_argument('--unique_id', type=str,
                        help='Unique identifier for the dataset. Used to name the output files.', default=generate_id())
    args = parser.parse_args()
    main(args)
