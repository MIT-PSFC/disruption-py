#!/usr/bin/env python3

import argparse
import json
import pickle

import joblib
import pandas as pd
from sklearn.metrics import ConfusionMatrixDisplay, fbeta_score

from disruption_py.ml.evaluate import DEFAULT_ORDER, THESIS_ORDER, eval_shots, predict
from disruption_py.ml.models import load_model
from disruption_py.core.utils.math import *

order_mapping = {
    "DEFAULT_ORDER": DEFAULT_ORDER,
    "THESIS_ORDER": THESIS_ORDER,
    "DATASET_ORDER": None,
}


def main(args):
    data = pd.read_csv(args.data_path)
    model = load_model(args.model_path)
    data = predict(model, data, order_mapping[args.order])
    data.to_csv(args.output_dir + f"eval_data_{args.unique_id}.csv")
    results = eval_shots(data, disruptivity=0.85)
    conf_mat = np.array(
        [[results["TN"], results["FP"]], [results["FN"], results["TP"]]]
    )
    print("Confusion Matrix Percentages", conf_mat / conf_mat.sum())
    print(
        "F2-Score:",
        5 * results["TP"] / (5 * results["TP"] + results["FP"] + 4 * results["FN"]),
    )
    disp = ConfusionMatrixDisplay(
        confusion_matrix=conf_mat, display_labels=model.classes_
    )
    disp.plot()
    save_open_plots(args.output_dir + f"eval_shots_{args.unique_id}.pdf")
    with open(args.output_dir + f"eval_{args.unique_id}.json", "w") as f:
        json.dump(vars(args), f)
    if args.visualize:
        if results["FP"] + results["TP"] + results["FN"] + results["TN"] < 20:
            plt.show()
        else:
            print(
                "Displaying plots would open too many figures (>20). Plots have been saved but will not be displayed"
            )
    print(f"Unique ID for this run: {args.unique_id}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # parser.add_argument('--model_path', type=str,
    #                     default='/fusion/projects/disruption_warning/papers/SORI2020/forest_245_20_dan.pkl')
    parser.add_argument(
        "--model_path",
        type=str,
        default="/fusion/projects/disruption_warning/papers/SORI2020/sori-repo/forest_245_15_dan.h5",
    )
    parser.add_argument("--output_dir", type=str, default="./output/")
    parser.add_argument("data_path", type=str)
    parser.add_argument("--visualize", type=bool, default=True)
    parser.add_argument(
        "--unique_id",
        type=str,
        help="Unique identifier for the dataset. Used to name the output files.",
        default=generate_id(),
    )
    parser.add_argument(
        "--order",
        type=str,
        default="THESIS_ORDER",
        choices=["DEFAULT_ORDER", "THESIS_ORDER", "DATASET_ORDER"],
    )
    args = parser.parse_args()
    main(args)
