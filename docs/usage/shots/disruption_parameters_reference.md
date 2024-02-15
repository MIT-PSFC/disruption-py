
## Disruption Parameter Descriptions

| Parameter | Description | Units | Validity Range |
|---|---|---|---|
| Greenwald_fraction | Greenwald density fraction = n_e/n_G, where n_e is the line-averaged density and n_G is Greenwald density limit | - | [0, 1.5] |
| Te_width | Half width half max of electron temperature profile from Thomson scattering | m | [0.04, 0.5] |
| Wmhd | Total magnetic energy stored in the plasma | J | [0, 2e5] |
| beta_p | Plasma poloidal beta, ratio between plasma pressure and magnetic pressure | - | [0, 1.1] |
| beta_n | Normalized beta, ratio between plasma kinetic energy and magnetic energy | - | [0, 2] |
| dipprog_dt | Time derivative of the programmed plasma current | A/s | - |
| intentional_disruption | Whether a disruption was unintentional (0), intentional (1), or non-disrupted (NaN) | - | [0,1,NaN] |
| ip | Plasma current | Ampere | - |
| ip_error | Error on the plasma current (ip-ipprog) | - | - |
| kappa | Plasma elongation | m | [0.8, 2] |
| li | Plasma normalized internal inductance | - | [0.2, 4.5] |
| lower_gap | Lower gap | m | [0.025, 0.3] |
| n_e | Line-averaged electron density of the plasma core | m^-3 | - |
| n_equal_1_mode | n=1 component of the perturbed magnetic field | Tesla | - |
| n_over_ncrit | Vertical stability parameter | - | [-0.5, 2] |
| p_icrf | Ion cyclotron power | W | [0, 6e6] |
| p_lh | Lower hybrid power | W | [0, 1e6] |
| p_oh | Ohmic power | W | [0, 20e6] |
| p_rad | Radiated power from the plasma | W | [0, 20e6] |
| q0 | Safety factor at the core plasma | - | [0, 10] |
| q95 | Safety factor at 95% of poloidal flux surface | - | [0, 20] |
| qstar | Cylindrical safety factor | - | [0, 30] |
| radiated_fraction | Total injected power from the beams divided by the radiated power from the plasma | - | [0,2-3] |
| shot | Discharge identifier, replicated per each time slice | - | - |
| ssep | Distance on midplane between 1st and 2nd separatrices | m | - |
| time | Time during the discharge | s | - |
| time_until_disrupt | Elapsing time before the disruption event. Target variable | s | NaN or numeric |
| upper_gap | Upper gap | m | [0, 0.21] |
| v_loop | Edge loop voltage; time derivative of a weighted average of flux loops obtained from MFLUXloop voltage | V | [-7, 26] |
| z_error | Difference between the actual position of the current centroid and the requested one (Z_prog) | m | - |
| zcur | Actual vertical position of the current centroid, z_error - Z_prog | m | - |
| Mirnov | Fluctuation amplitude of one magnetic probe. Measurement of MHD activity and plasma instability | Tesla/s | [0, 50] |