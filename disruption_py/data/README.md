
# C-MOD and DIII-D shotlists

### CMod_vde.csv

This is a list of disruptive shots with Vertical Displacement Events (VDEs) on Alcator C-Mod.

These shots were taken between 2012 -> 2016.

This list does not descriminate between hot and cold VDEs. Where "hot" refers to the VDEs preceding the Thermal Quench (TQ), while "old" refers to the VDEs following it.

This file contains three columns, "shot", "onset_time_seconds", and "notes". The "shot" column contains the shot numbers, the "onset_time_seconds" contains the manually labeled onset times of the VDEes in units of seconds, and the "notes" column contains additional notes about the shot (if present). Shot numbers on C-Mod also correspond to the date and order in which the shot was taken, following the pattern (1)YY-MM-DD-XXX, where YY is the last 2 digits of the year (with "1" being included for 2000+) and XXX is the XXXth shot taken that day.

All of these shots have been manually vetted (by @AlexSaperstein) for the presence of VDEs and their onset times. Onset times were chosen based on significant deviations (> ~2 cm) in the "z_error" feature.

It should be noted that onset times were chosen using a version of the disruption-py dataset that had poor time-resolution more than 20ms before the disruption. As a result, onset times which precede the disruption by more than 20ms have an error of ~20ms associated with them.