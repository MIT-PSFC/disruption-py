
## C-Mod Disruption Parameter Descriptions

| Parameter | Description | Units | Validity Range |
|---|---|---|---|
| a_minor | Minor radius (EFIT). | m | - |
| beta_n | Normalized plasma beta computed from EFIT. | dimensionless | [0, 2] |
| beta_p | Beta poloidal (EFIT). | dimensionless | [0, 1.5] |
| bt | Toroidal magnetic field strength. | T | - |
| chisq | Magnetic chi^2 (EFIT). | dimensionless | [0.1, 40] |
| dbetap_dt | Time derivative of the poloidal beta (EFIT). | s^-1 | - |
| dip_dt | Time derivative of the actual plasma current. | A/s | - |
| dip_smoothed | Time derivative of the smoothed measured plasma current. | A/s | - |
| dipprog_dt | Time derivative of the programmed plasma current. | A/s | [-5e6, 5e6] |
| dli_dt | Time derivative of the internal inductance (EFIT). | s^-1 | - |
| dn_dt | Time derivative of the measured line-averaged electron density. | m^-3*s^-1 | - |
| dprad_dt | Time derivative of the measured total radiated power. | W/s | - |
| dwmhd_dt | Time derivative of the total plasma energy (EFIT). | J/s | - |
| greenwald_fraction | Greenwald fraction. | dimensionless | [0, 1.5] |
| h98 | H98 energy confinement time enhancement factor over the IPB98(y,2) scaling. | dimensionless | [0.05, 1.75] |
| h_alpha | H-alpha line emission intensity. | W/(m^2*sr) | [0, 161.921] |
| i_efc | Error field correction (EFC) coil current. | A | - |
| ip | Measured plasma current (Rogowski coil). | A | - |
| ip_error | Difference between the measured and programmed plasma current. | A | - |
| ip_prog | Programmed plasma current. | A | - |
| kappa | Elongation at plasma boundary (EFIT). | dimensionless | [0.8, 2] |
| kappa_area | Elongation calculated using the plasma area and minor radius from EFIT. | dimensionless | [0.8, 2] |
| lh_power_threshold | Martin 2008 L-H transition power threshold scaling. | W | [0, 1e7] |
| li | Normalized internal inductance, li(3) (EFIT). | dimensionless | [0.2, 4.5] |
| lower_gap | Bottom gap between x-point and wall (EFIT). | m | [0.025, 0.25] |
| n_e | Measured line-averaged electron density. | m^-3 | [2e19, 7e20] |
| n_equal_1_mode | Amplitude of the n=1 mode. | T | [0, 0.02] |
| n_equal_1_normalized | Amplitude of the n=1 mode normalized to the toroidal magnetic field. | dimensionless | - |
| n_equal_1_phase | Phase of the n=1 mode. | rad | - |
| n_over_ncrit | Vertical stability parameter, vacuum field index normalized to critical index. (EFIT). | dimensionless | [-3, 3] |
| ne_peaking | Peaking factor of the electron density profile measured by the Thomson scattering diagnostic. | dimensionless | - |
| p_icrf | Total ion cyclotron resonant heating power. | W | [0, inf] |
| p_input | Total input power (ohmic + LH + ICRF). | W | - |
| p_lh | Total lower hybrid heating power. | W | [0, inf] |
| p_oh | Total Ohmic heating power. | W | [0, 5e6] |
| p_rad | Total radiated power measured by the 2pi diode. | W | [0, inf] |
| prad_peaking | Peaking factor of the radiated power. | dimensionless | - |
| pressure_peaking | Peaking factor of the pressure profile measured by the Thomson scattering diagnostic. | dimensionless | - |
| q0 | Safety factor at magnetic axis (EFIT). | dimensionless | [0, 2] |
| q95 | Safety factor at 95% flux surface (EFIT). | dimensionless | [0, inf] |
| qstar | Equivalent (kink) safety factor (EFIT). | dimensionless | [0, inf] |
| radiated_fraction | Total radiated power fraction. | dimensionless | [0, 15] |
| rmagx | Major radius of magnetic axis (EFIT). | m | - |
| ssep | Outboard radial distance to external second separatrix for single null configurations. Positive for single top-null and negative for single bottom-null. Near-zero for double-null (EFIT). | m | [-1, 1] |
| sxr | Core soft X-ray measurement. | W | - |
| tau_rad | Radiative cooling time-scale. | s | - |
| te_core_vs_avg_ece | Core-vs-average peaking factor of the electron temperature profile from the electron cyclotron emission (ECE) data. | dimensionless | - |
| te_edge_vs_avg_ece | Edge-vs-average peaking factor of the electron temperature profile from the electron cyclotron emission (ECE) data. | dimensionless | - |
| te_peaking | Peaking factor of the electron temperature profile measured by the Thomson scattering diagnostic. | dimensionless | - |
| te_width | Half-width at half-maximum of the Gaussian fit to the core electron temperature profile. | m | - |
| te_width_ece | Half-width at half-maximum of the Gaussian fit to the electron temperature profile from the electron cyclotron emission (ECE) data. | m | - |
| thermal_quench_time | Time of the disruptive thermal quench onset. NaN for non-disruptive discharges. | s | [0, 3] |
| time_until_disrupt | Time until the current quench. NaN for non-disruptive discharges. | s | - |
| tribot | Bottom triangularity (EFIT). | dimensionless | - |
| tritop | Top triangularity (EFIT). | dimensionless | - |
| upper_gap | Top gap between x-point and wall (EFIT). | m | [0, inf] |
| v_loop | Measured loop voltage. | V | [-7, 7] |
| v_loop_efit | Reconstructed loop voltage (EFIT), smoothed using non-causal filter. | V | [-7, 7] |
| v_surf | Plasma surface voltage computed from EFIT. | V | - |
| v_z | Vertical velocity of the plasma current centroid. | m/s | - |
| wmhd | Total plasma energy (EFIT). | J | [5000, 300000] |
| z_error | Difference between the measured and programmed vertical position of the plasma current centroid. | m | [-1, 1] |
| z_prog | Programmed vertical position of the plasma current centroid. | m | - |
| z_times_v_z | Product of the vertical position and vertical velocity of the plasma current centroid. | m^2/s | - |
| zcur | Measured vertical position of the plasma current centroid. | m | [-0.5, 0.5] |

## DIII-D Disruption Parameter Descriptions

| Parameter | Description | Units | Validity Range |
|---|---|---|---|
| aminor | Minor radius (EFIT). | m | - |
| beta_n | Normalized plasma beta computed from EFIT. | dimensionless | [0, 2] |
| beta_p | Beta poloidal (EFIT). | dimensionless | [0, 1.5] |
| beta_p_rt | Beta poloidal from real-time EFIT calculation. | dimensionless | [0, 1.5] |
| btor | Toroidal magnetic field strength. | T | [-2.5, 2.5] |
| chisq | Magnetic chi^2 (EFIT). | dimensionless | [0.1, 70] |
| chisq_rt | Magnetic chi^2 from real-time EFIT calculation. | dimensionless | [0.1, 70] |
| current_quench_time | Time of the current quench event. NaN for non-disruptive discharges. | s | - |
| dbetap_dt | Time derivative of the poloidal beta (EFIT). | s^-1 | - |
| dbetap_dt_rt | Time derivative of the poloidal beta from real-time EFIT calculation. | s^-1 | - |
| delta | Triangularity at plasma boundary (EFIT). | dimensionless | - |
| dip_dt | Time derivative of the actual plasma current. | A/s | - |
| dip_dt_rt | Time derivative of the real-time actual plasma current. | A/s | - |
| dipprog_dt | Time derivative of the programmed plasma current. | A/s | [-5e6, 5e6] |
| dipprog_dt_rt | Time derivative of the real-time programmed plasma current. | A/s | [-5e6, 5e6] |
| dli_dt | Time derivative of the internal inductance (EFIT). | s^-1 | - |
| dn_dt | Time derivative of the measured line-averaged electron density. | m^-3*s^-1 | - |
| dn_dt_rt | Time derivative of the measured line-averaged electron density available to PCS in real-time. | m^-3*s^-1 | - |
| dwmhd_dt | Time derivative of the total plasma energy (EFIT). | J/s | - |
| greenwald_fraction | Greenwald fraction. | dimensionless | [0, 1.5] |
| greenwald_fraction_rt | Greenwald fraction calculated from real-time density. | dimensionless | [0, 1.5] |
| h98 | H98 energy confinement time enhancement factor over the IPB98(y,2) scaling. | dimensionless | [0.05, 1.75] |
| h_alpha | H-alpha line emission intensity. | W/(m^2*sr) | [0, 161.921] |
| ip | Measured plasma current (Rogowski coil). | A | - |
| ip_error | Difference between the measured and programmed plasma current. | A | - |
| ip_error_rt | Difference between the real-time measured and programmed plasma current. | A | - |
| ip_prog | Programmed plasma current. | A | - |
| ip_prog_rt | Real-time programmed plasma current. | A | - |
| ip_rt | Real-time measured plasma current (Rogowski coil). | A | - |
| kappa | Elongation at plasma boundary (EFIT). | dimensionless | [0.8, 2] |
| kappa_area | Elongation calculated using the plasma area and minor radius from EFIT. | dimensionless | [0.8, 2] |
| li | Normalized internal inductance, li(3) (EFIT). | dimensionless | [0.2, 4.5] |
| li_rt | Normalized internal inductance, li(3) from real-time EFIT calculation. | dimensionless | [0.2, 4.5] |
| lower_gap | Bottom gap between x-point and wall (EFIT). | m | [0.025, 0.25] |
| n1rms | RMS amplitude of the n=1 magnetic field perturbation. | T | [0, 0.02] |
| n1rms_normalized | RMS of the n=1 mode normalized to the toroidal magnetic field strength. | dimensionless | - |
| n_e | Measured line-averaged electron density. | m^-3 | [2e19, 7e20] |
| n_e_rt | Real-time measured line-averaged electron density available to PCS. | m^-3 | [2e19, 7e20] |
| n_equal_1_mode | Amplitude of the n=1 mode. | T | [0, 0.02] |
| n_equal_1_normalized | Amplitude of the n=1 mode normalized to the toroidal magnetic field. | dimensionless | - |
| ne_peaking_cva_rt | Peaking factor of the electron density profile, core vs all channels. | dimensionless | - |
| p_ech | Total electron cyclotron resonant heating power. | W | [0, inf] |
| p_nbi | Total neutral beam injection heating power. | W | [0, inf] |
| p_ohm | Total Ohmic heating power. | W | [0, 5e6] |
| p_rad | Total radiated power measured by 2pi foil bolometer. | W | [0, inf] |
| power_supply_railed | Power supply railed signal indicator for PCS feedback control of Ip. | dimensionless | [0, 1] |
| prad_peaking_cva_rt | Peaking factor for the radiated power, core vs all-but-divertor channels. | dimensionless | - |
| prad_peaking_xdiv_rt | Peaking factor for the radiated power, divertor vs all-but-core channels. | dimensionless | - |
| q0 | Safety factor at magnetic axis (EFIT). | dimensionless | [0, 2] |
| q95 | Safety factor at 95% flux surface (EFIT). | dimensionless | [0, inf] |
| q95_rt | Safety factor at 95% flux surface from real-time EFIT calculation. | dimensionless | [0, inf] |
| qstar | Equivalent (kink) safety factor (EFIT). | dimensionless | [0, inf] |
| radiated_fraction | Total radiated power fraction. | dimensionless | [0, 15] |
| squareness | Average of lower-outer and upper-outer squareness, T. Luce PPCF 55 (2013) 095009. | dimensionless | - |
| te_peaking_cva_rt | Peaking factor of the electron temperature profile, core vs all channels. | dimensionless | [0, inf] |
| time_until_disrupt | Time until the current quench. NaN for non-disruptive discharges. | s | - |
| upper_gap | Top gap between x-point and wall (EFIT). | m | [0, inf] |
| v_loop | Measured loop voltage. | V | [-7, 7] |
| wmhd | Total plasma energy (EFIT). | J | [5000, 300000] |
| wmhd_rt | Total plasma energy from real-time EFIT calculation. | J | [5000, 300000] |
| z_eff | Effective charge. | dimensionless | - |
| zcur | Measured vertical position of the plasma current centroid. | m | [-0.5, 0.5] |
| zcur_normalized | Measured vertical position of the plasma current centroid normalized to plasma minor radius. | dimensionless | - |

## Generic Disruption Parameter Descriptions

| Parameter | Description | Units | Validity Range |
|---|---|---|---|
| current_quench_time | Time of the current quench event. NaN for non-disruptive discharges. | s | - |
| time_domain | Categorical phase of shot at each time. 1: ramp-up, 2: flat-top, 3: ramp-down. | dimensionless | - |

## MAST Disruption Parameter Descriptions

| Parameter | Description | Units | Validity Range |
|---|---|---|---|
| a_minor | Minor radius of the plasma boundary (defined as (Rmax-Rmin) / 2 of the boundary) | m | - |
| beta_n | Toroidal beta, defined as the volume-averaged total perpendicular. | dimensionless | - |
| beta_p | Beta poloidal (EFIT). | dimensionless | - |
| beta_t | Beta toroidal (EFIT). | dimensionless | - |
| bphi_rmag | Bphi at rmag | T | - |
| bvac_rmag | Bvac at rmag | T | - |
| d_alpha | D-alpha camera signal | V | - |
| dip_dt | Time derivative of the Actual plasma current. | A/s | - |
| dipprog_dt | Time derivative of the programmed plasma current. | A/s | - |
| dn_dt | Time derivative of the measured line-averaged electron density. | m^-3*s^-1 | - |
| gas_inboard_total | Gas injected by the inboard valves | cps | - |
| gas_outboard_total | Gas injected by the outboard valves | cps | - |
| gas_total_injected | Total cumulative injected gas | counts | - |
| greenwald_fraction | Greenwald fraction. | dimensionless | - |
| ip | Actual plasma current. | A | - |
| ip_prog | Programmed plasma current. | A | - |
| kappa | Elongation at plasma boundary (EFIT). | dimensionless | - |
| li | Internal inductance (EFIT). | dimensionless | - |
| n_e | Measured line-averaged electron density. | m^-3 | - |
| ne_core | Core electron density measurement | m^-3 | - |
| ne_peaking | Peaking factor of the electron density profile measured by the Thomson scattering diagnostic. | dimensionless | - |
| p_nbi | NBI heating power. | W | - |
| p_oh | Total Ohmic heating power. | W | - |
| p_rad | Radiated power measured by the 2pi diode. | W | - |
| prad_peaking | Peaking factor of the radiated power. | dimensionless | - |
| pressure_peaking | Peaking factor of the electron pressure profile measured by the Thomson scattering diagnostic. | dimensionless | - |
| q95 | q at 95% flux surface (EFIT). | dimensionless | - |
| rmagx | Major radius of the magnetic axis | m | - |
| rmagz | Height of the magnetic axis | m | - |
| sxr | Core soft X-ray measurement. | V | - |
| sxr_core | Core soft X-ray measurement. | V | - |
| sxr_edge | Edge soft X-ray measurement. | V | - |
| te_core | Core electron temperature measurement | eV | - |
| te_peaking | Peaking factor of the electron temperature profile measured by the Thomson scattering diagnostic. | dimensionless | - |
| te_width | Half-width at half-maximum of the Gaussian fit to the electron temperature profile. | m | - |
| tribot | Bottom triangularity (EFIT). | dimensionless | - |
| tritop | Top triangularity (EFIT). | dimensionless | - |
| v_loop_dynamic | Dynamic LCFS loop voltage from EFIT. Defined as V_loop(t) = loop int(B_p * V_point(t)) / loop int(B_p)(t) where V_point = -2 * pi * d(PSI at each point) / dt | V | - |
| v_loop_static | Static LCFS loop voltage from EFIT. Defined as V_loop = -2 * pi * d(psi at LCFS) / dt | V | - |
| v_z | Vertical velocity of the plasma current centroid. | m/s | - |
| volume | Volume of the plasma. | m^3 | - |
| wmhd | Total plasma energy (EFIT). | J | - |
| z_error | Difference between the measured and programmed vertical position of the plasma current centroid. | m | - |
| z_prog | Programmed vertical position of the plasma current centroid. | m | - |
| z_times_v_z | Product of the vertical position and vertical velocity of the plasma current centroid. | m^2/s | - |
| zcur | Measured vertical position of the plasma current centroid. | m | - |
