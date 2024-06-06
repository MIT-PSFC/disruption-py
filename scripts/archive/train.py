#!/usr/bin/env python3

import argparse
import json
import pickle
from datetime import date

import joblib
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import ConfusionMatrixDisplay, confusion_matrix, fbeta_score

from disruption_py.ml.train import grid_search, train_local
from disruption_py.utils.math_utils import generate_id, save_open_plots


def main(args):
    if "pkl" in args.train_path:
        with open(args.train_path, "rb") as f:
            train = pickle.load(f, encoding="latin1")
        x_train = train["X_train"]
        y_train = train["y_train"].astype(int)
    else:
        train = pd.read_csv(args.train_path)
        y_train = train["label"]
        x_train = train.drop(
            ["label", "time_until_disrupt", "shot"], axis=1, errors="ignore"
        )
    if "pkl" in args.test_path:
        with open(args.test_path, "rb") as f:
            test = pickle.load(f, encoding="latin1")
        x_test = test["X_test"]
        y_test = test["y_test"]
    else:
        test = pd.read_csv(args.test_path)
        y_test = test["label"]
        x_test = test.drop(
            ["label", "time_until_disrupt", "shot"], axis=1, errors="ignore"
        )
    if args.grid_search is not None:
        grid = json.loads(args.grid_search)
        results = grid_search(x_train, y_train, x_test, y_test, **grid)
        with open(args.output_dir + f"grid_search_{args.unique_id}_.json", "w") as f:
            json.dump(results, f)
        return
    else:
        results = [
            train_local(
                x_train,
                y_train,
                x_test,
                y_test,
                features=args.features,
                model=args.model,
            )
        ]
        best_result = results[0]
    model = best_result["model"]
    train["scores"] = best_result["predictions_proba"][:, 1]
    test["scores"] = model.predict(x_test)
    print("Train F2-Score:", fbeta_score(y_train, best_result["predictions"], beta=2))
    print("Test F2-Score:", fbeta_score(y_test, test["scores"], beta=2))
    conf_mat = confusion_matrix(y_train, best_result["predictions"])
    disp_train = ConfusionMatrixDisplay(
        confusion_matrix=conf_mat, display_labels=model.classes_
    )
    conf_mat = confusion_matrix(y_test, test["scores"])
    disp_test = ConfusionMatrixDisplay(
        confusion_matrix=conf_mat, display_labels=model.classes_
    )
    disp_train.plot()
    disp_test.plot()
    disp_train.ax_.set_title("Train Matrix")
    disp_test.ax_.set_title("Test Matrix")
    # Plot feature importances
    plt.figure()
    if args.model == "LogisticRegression":
        feature_names = model.named_steps["polynomialfeatures"].get_feature_names(
            input_features=x_train.columns
        )
    else:
        feature_names = x_train.columns
    feature_importances = zip(feature_names, best_result["features"])
    sorted_features = sorted(feature_importances, key=lambda x: x[1], reverse=True)
    top_features = sorted_features[:10]
    top_feature_names, top_importances = zip(*top_features)
    plt.barh(top_feature_names, top_importances)
    plt.title("Feature Importances")
    plt.xlabel("Relative Importance")
    plt.ylabel("Feature")
    plt.show()
    # Save plots using ../disruption_py/utils
    save_open_plots(f"train_plots_{args.unique_id}.pdf")
    joblib.dump(
        model, args.output_dir + f"{args.model}{date.today()}_{args.unique_id}.joblib"
    )
    with open(args.output_dir + f"train_{args.unique_id}_.json", "w") as f:
        json.dump(vars(args), f)
    print(f"Unique ID for this run: {args.unique_id}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train a model on a dataset.")
    parser.add_argument(
        "--train_path",
        type=str,
        help="Path to training data. Must be a csv file.",
        default="./output/train.csv",
    )
    parser.add_argument(
        "--test_path",
        type=str,
        help="Path to test data. Must be a csv file.",
        default="./output/test.csv",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="./output/",
        help="Path to output model files. Must be a directory.",
    )
    parser.add_argument(
        "--grid_search",
        type=str,
        default=None,
        help="Path to grid search parameters. Must be a json file.",
    )
    parser.add_argument(
        "--features",
        type=str,
        nargs="+",
        help="List of features to use for training. If not provided, all features will be used.",
        default=None,
    )
    parser.add_argument(
        "--unique_id",
        type=str,
        help="Unique identifier for the dataset. Used to name the output files.",
        default=generate_id(),
    )
    parser.add_argument(
        "--model",
        type=str,
        default="RandomForestClassifier",
        choices=[
            "RandomForestClassifier",
            "LogisticRegression",
            "SupportVectorMachine",
            "LinearClassifier",
        ],
        help="Model to use for training. Must be a valid sklearn model.",
    )
    args = parser.parse_args()
    main(args)
