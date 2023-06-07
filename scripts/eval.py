import argparse
import pickle
import joblib
import json

import pandas as pd
from sklearn.metrics import ConfusionMatrixDisplay

from disruption_py.utils import *
from disruption_py.ml.models import load_model 
from disruption_py.ml.evaluate import eval_shots, predict, DEFAULT_ORDER

def main(args):
    data = pd.read_csv(args.data_path)
    model = load_model(args.model_path)
    data = predict(model, data)
    data.to_csv(args.output_dir + f"eval_data_{args.unique_id}.csv")
    results = eval_shots(data)
    print(results)
    conf_mat = np.array([[results['TN'], results['FP']],
                         [results['FN'], results['TP']]])
    disp = ConfusionMatrixDisplay(
        confusion_matrix=conf_mat, display_labels=model.classes_)
    disp.plot()
    save_open_plots(args.output_dir + f"eval_shots_{args.unique_id}.pdf")
    with open(args.output_dir + f"eval_{args.unique_id}.json", "w") as f:
        json.dump(vars(args), f)
    if args.visualize:
        if results['FP'] + results['TP'] + results['FN'] + results['TN'] < 20:
            plt.show()
        else:
            print("Displaying plots would open too many figures (>20). Plots have been saved but will not be displayed")
    print(f"Unique ID for this run: {args.unique_id}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--model_path', type=str,
                        default='/fusion/projects/disruption_warning/papers/SORI2020/forest_245_20_dan.pkl')
    parser.add_argument('--output_dir', type=str, default='./output/')
    parser.add_argument('data_path', type=str)
    parser.add_argument('--visualize', type=bool, default=True)
    parser.add_argument('--unique_id', type=str,
                        help='Unique identifier for the dataset. Used to name the output files.', default=generate_id())
    args = parser.parse_args()
    main(args)
