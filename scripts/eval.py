import argparse
import pickle
import joblib
import json

import pandas as pd
from sklearn.metrics import confusion_matrix, fbeta_score, ConfusionMatrixDisplay
from sklearn.ensemble import RandomForestClassifier
from disruption_py.utils import *

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
                print(f"Warning time: {warning_time}")
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
    df['score'] = model.predict_proba(arr)
    return df


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
        plt.show()
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
