from disruption_py.utils import *


def eval_shots(df, lower_threshold = .05, disruptivity = .45, window= .025):
    for _, shot_data in df.groupby(by=['shot']):
        alarm = trigger_alarm(shot_data['time'], shot_data['score'],lower_threshold,disruptivity, window)
   
