
# C-MOD shotlists

> [!TIP]
> Shot numbers on Alcator C-Mod correspond to the date and order in which the shot was taken, with the following the pattern: `[1]YYMMDDNNN`, where `YY` is the last 2 digits of the year (`99` for 1999, `100` for 2000, and so on), `MM` is the month, `DD` is the day, and `NNN` is the N-th shot taken that day.

## Vertical Displacement Events (VDEs)

info|value
-|-
_file_ | `cmod_vde.csv`
_machine_ | **Alcator C-MOD**
_total shots_ | **99**
_years_ | **2012-2016**
_maintainer_ | @AlexSaperstein

This is a list of disruptive shots with Vertical Displacement Events (VDEs) on Alcator C-Mod.

This list does not descriminate between hot and cold VDEs.
Where 'hot' refers to the VDEs preceding the Thermal Quench (TQ), while 'cold' refers to the VDEs following it.

### Details

Column description:
- `shot` (integer): shot numbers,
- `onset_time` (seconds): manually labeled onset times of the VDEs,
- `notes` (string): any additional notes about the shot.

All of these shots have been manually vetted by Alex Saperstein for the presence of VDEs and their onset times.
Onset times were chosen based on significant deviations (> ~2 cm) in the `z_error` feature.

> [!NOTE]
> It should be noted that onset times were chosen using a version of the disruption-py dataset that had poor time-resolution more than 20ms before the disruption.
> As a result, onset times which precede the disruption by more than 20ms have an error of ~20ms associated with them.

## "Unidentified Flying Object" (UFO) disruptions

info|value
-|-
_file_ | `cmod_ufo.csv`
_machine_ | **Alcator C-MOD**
_total shots_ | **122**
_years_ | **2012-2016**
_maintainer_ | @HWietfeldt

This is a list of disruptive shots with likely UFO disruptions on Alcator C-Mod.

This list was generated from an automatically generated database of high-Z injection events on C-Mod discharges.
High-Z injection events were labeled using a Mo +31 charge state signal measured from a VUV spectrometer.
Increases in the Mo +31 charge state signal were labeled as injections if:
- global radiated power increased by 1 MW around the time of injection and before the thermal quench,
- soft x-ray radiation increased around the time of injection and before the thermal quench,
- low error in the vertical position at the time of injection.

This list contains likely UFO disruptions, defined as shots with a labeled high-Z injection within 10 ms of the thermal quench.
Only 1/6th of these shots were manually validated as containing a UFO disruption by Henry Wietfeldt.
Whether a shot has been validated is indicated in the "notes" column.
Other shots may have labeling errors since they have not yet been manually validated.

### Details

Column description:
- `shot` (integer): shot numbers,
- `injection_time` (seconds): time of the injection, based on when the Mo +31 charge state signal began to increase,
- `duration` (seconds): time between when the Mo +31 charge state signal began to increase and stopped increasing,
- `t_until_thermal_quench` (seconds): time between `t_injection` and the onset of the thermal quench, determined from core Te measurements using an ECE diagnostic system,
- `notes` (string): additional notes, including whether the shot has been validated as a UFO.
