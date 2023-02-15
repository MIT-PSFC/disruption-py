import argparse
import pickle
import joblib

import pandas as pd
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, fbeta_score

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
    'time': 8
}


def eval_shots(df, lower_threshold=.05, disruptivity=.45, window=.025):
    good_warnings = []
    missed_warnings = []
    false_alarms = []
    warning_times = dict()
    disruptions = 0
    non_disruptions = 0
    for _, shot_data in df.groupby(by=['shot']):
        shot_id = shot_data['shot'].values[0]
        time_until_disrupt = shot_data['time_until_disrupt'].values
        alarm = trigger_alarm(
            shot_data['time'], shot_data['score'], lower_threshold, disruptivity, window)
        trigger = np.where(np.array(alarm) == 1)[0]
        disrupted = np.isnan(time_until_disrupt).all()
        disruptions += disrupted
        non_disruptions += not disrupted
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
    print(
        "Disruptions warned: %d/%d (%d%%)"
        % (len(good_warnings), len(disruptions), round(float(len(good_warnings)) / len(disruptions) * 100, 2))
    )
    print(
        "Missed disruptions: %d/%d (%d%%)"
        % (len(missed_warnings), len(disruptions), round(float(len(missed_warnings)) / len(disruptions) * 100, 2))
    )
    print(
        "False Alarms: %d/%d (%d%%)"
        % (len(false_alarms), len(non_disruptions), round(float(len(false_alarms)) / len(non_disruptions) * 100, 2))
    )

    results['good_warnings'] = good_warnings
    results['missed_warnings'] = missed_warnings
    results['false_alarms'] = false_alarms
    results['warning_times'] = warning_times
    results['TP'] = len(good_warnings) / float(len(disruptions))
    results['FP'] = len(false_alarms) / float(len(non_disruptions))
    results['FN'] = len(missed_warnings) / float(len(disruptions))
    results['TN'] = -(len(good_warnings) + len(false_alarms) +
                      len(missed_warnings) - len(df) / float(len(non_disruptions)))
    return results


def load_model(model_path):
    if model_path.endswith('.pkl'):
        model = pickle.load(open(model_path, 'rb'))
    elif model_path.endswith('.joblib'):
        model = joblib.load(model_path)
    elif model_path.endswith('.h5'):
        pass
    else:
        raise ValueError(f"Unknown model file type: {model_path}")
    return model


def predict(model, df, order=DEFAULT_ORDER):
    arr = np.empty((len(df), len(order)))
    for col in order:
        arr[:, order[col]] = df[col].values
    df['score'] = model.predict(arr)
    return df


def main(args):
    data = pd.read_csv(args.data_path)
    model = load_model(args.model_path)
    data = predict(model, data)
    results = eval_shots(data)
    conf_mat = np.array([[results['TN'], results['FP']],
                        [results['FN'], results['TP']]])
    disp_train = ConfusionMatrixDisplay(
        confusion_matrix=conf_mat, display_labels=model.classes_)
    disp_train.plot()
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('model_path', type=str)
    parser.add_argument('data_path', type=str)
