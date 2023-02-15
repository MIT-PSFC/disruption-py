import argparse
import pickle
import joblib

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


def plot_confusion_matrix(cm, classes=[0, 1],
                          normalize=False,
                          title=None,
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    print(cm)

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           xticklabels=classes, yticklabels=classes,
           title=title,
           ylabel='True label',
           xlabel='Predicted label')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    return ax


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
            shot_data['time'], shot_data['score'], lower_threshold, disruptivity, window)
        plt.figure()
        plt.plot(shot_data['time'], shot_data['score'])
        plt.plot(shot_data['time'], alarm)
        trigger = np.where(np.array(alarm) == 1)[0]
        disrupted = np.isnan(time_until_disrupt).all()
        if disrupted:
            disruptions.append(shot_id)
        else:
            non_disruptions.append(shot_id)
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
    if len(disruptions) == 0:
        print(
            "Disruptions warned: %d/%d (%d%%)"
            % (len(good_warnings), len(disruptions), round(float(len(good_warnings)) / len(disruptions) * 100, 2))
        )
        print(
            "Missed disruptions: %d/%d (%d%%)"
            % (len(missed_warnings), len(disruptions), round(float(len(missed_warnings)) / len(disruptions) * 100, 2))
        )
    if len(non_disruptions) == 0:
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
    save_open_plots('eval_shots.pdf')
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
    results = eval_shots(data)
    print(results)
    conf_mat = np.array([[results['TN'], results['FP']],
                         [results['FN'], results['TP']]])
    disp = ConfusionMatrixDisplay(
        confusion_matrix=conf_mat, display_labels=model.classes_)
    disp.plot()
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--model_path', type=str,
                        default='/fusion/projects/disruption_warning/papers/SORI2020/forest_245_20_dan.pkl')
    parser.add_argument('data_path', type=str)
    args = parser.parse_args()
    main(args)
