import numpy as np

from disruption_py.utils import *
from disruption_py.ml.models import load_model_from_hdf5

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

def predict(model, df, order=DEFAULT_ORDER):
    arr = np.empty((len(df), len(order)))
    for col in order:
        arr[:, order[col]] = df[col].values
    predictions = model.predict_proba(arr)
    if len(predictions.shape) > 1:
        predictions = predictions[:,1]
    df['score'] = predictions
    return df

def hist_time_th(thL, thH, thA, dt, v, v_old, c_old, a_old):

    # =================
    # double threshold
    # =================
    # if signal is below LOW threshold output is always 0
    if v < thL:
        v = 0
    # if signal is above HIGHT threshold output is always 1
    elif v > thH:
        v = 1
    else:
        v = v_old
    # =================
    # time accumulation
    # =================
    # when output signal is 0 counter gets reset
    if v == 0:
        c = 0
    # when output signal is 1 counter gets incremented
    else:
        c = c_old + dt
    # =================
    # sound alarm
    # =================
    # if counter surpasses allarm threshold or old allarm was on
    if c > thA or a_old:
        a = 1
    else:
        a = 0

    return v, c, a


def trigger_alarm(time, value, low_thr, high_thr, thA):
    output = []
    alarm = []
    if time[0] < 0:
        tstart = 0
    else:
        tstart = time[0]  # 0
    v_old = 0
    c_old = 0
    a_old = 0

    for t, v in zip(time, value):
        dt = t - tstart
        v, c, a = hist_time_th(low_thr, high_thr, thA,
                               dt, v, v_old, c_old, a_old)
        output.append(v)
        alarm.append(a)
        if output[0] != 0:
            v, c, a = 0, 0, 0
            output[0], alarm[0] = 0, 0
        v_old, c_old, a_old = v, c, a
        tstart = t
    return alarm
