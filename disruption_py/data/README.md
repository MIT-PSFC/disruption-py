
# C-MOD and DIII-D shotlists

### cmod_vde.csv

Maintainer: @AlexSaperstein

This is a list of disruptive shots with Vertical Displacement Events (VDEs) on Alcator C-Mod.\
These shots were taken between 2012 -> 2016.\
This list does not descriminate between hot and cold VDEs. Where "hot" refers to the VDEs preceding the Thermal Quench (TQ), while "old" refers to the VDEs following it.
This file contains three columns, "shot", "onset_time_seconds", and "notes". The "shot" column contains the shot numbers, the "onset_time_seconds" contains the manually labeled onset times of the VDEes in units of seconds, and the "notes" column contains additional notes about the shot (if present). Shot numbers on C-Mod also correspond to the date and order in which the shot was taken, following the pattern (1)YY-MM-DD-XXX, where YY is the last 2 digits of the year (with "1" being included for 2000+) and XXX is the XXXth shot taken that day.\
All of these shots have been manually vetted (by @AlexSaperstein) for the presence of VDEs and their onset times. Onset times were chosen based on significant deviations (> ~2 cm) in the "z_error" feature.

> [!NOTE]\
> It should be noted that onset times were chosen using a version of the disruption-py dataset that had poor time-resolution more than 20ms before the disruption. As a result, onset times which precede the disruption by more than 20ms have an error of ~20ms associated with them.

### CMod_ufo.csv

This is a list of disruptive shots with likely UFO disruptions on Alcator C-Mod.

These shots were taken between 2012->2016.

This list was generated from an automatically generated database of high-Z injection events on C-Mod discharges.
High-Z injection events were labeled using a Mo +31 charge state signal measured from a VUV spectrometer.
Increases in the Mo +31 charge state signal were labeled as injections if
- Global radiated power increased by 1 MW around the time of injection and before the thermal quench
- Soft x-ray radiation increased around the time of injection and before the thermal quench
- Low error in the vertical position at the time of injection

This list contains likely UFO disruptions, defined as shots with a labeled high-Z injection within 10 ms of the thermal quench. Only 1/6 of these shots were manually validated @hwietfeldt as containing a UFO disruption by. Whether a shot has been validated is indicated in the "notes" column. Other shots may have labeling errors since they have not yet been manually validated.

This file contains five columns:
- "shot": The shot number, following the pattern (1)YY-MM-DD-XXX, where YY is the last 2 digits of the year (with "1" being included for 2000+) and XXX is the XXXth shot taken that day
- "t_injection_seconds": The time (in seconds) of the injection, based on when the Mo +31 charge state signal began to increase
- "duration_seconds": The time (in seconds) between when the Mo +31 charge state signal began to increase and stopped increasing.
- "t_until_thermal_quench_seconds": The time (in seconds) between "t_injection_seconds" and the onset of the thermal quench, determined from core Te measurements using an ECE diagnostic system
- "notes": Additional notes, including whether the shot has been validated as a UFO.

For questions, contact Henry Wietfeldt (Github: @hwietfeldt, Email: henrycw@mit.edu)

