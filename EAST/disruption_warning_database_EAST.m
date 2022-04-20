function disruption_warning_database_EAST(shotlist);

% 'shotlist' is a Matlab array of shot numbers.  It can be a single shot,
%   or a sequence of shots [49848:51390], or a list of specific shots
%   [52708, 52709, 53277], or combinations of these [49848:51390, 52708,
%   52709, 53277].
%
% Examples of usage:
%
%  disruption_warning_database(57508);  % Process shot 57508 only
%
%  disruption_warning_database(49848:51390);  % Process all plasma shots
%                                             % sequentially from 49848 
%                                             % to 51390.
% Revision/history:
%  2015/12    - RSG; Created main routine (adapted from the disruption
%                     warning database program that I wrote for the C-Mod
%                     tokamak) 
%  2015/12/10 - RSG: Implemented the routine "get_Ip_parameters"
%  2015/12/11 - RSG; Wrote code segment to implement non-uniform time
%                     sampling, with much more frequent sampling just prior
%                     to disruption times.
%  2015/12/15 - RSG; Implemented the routines "get_Z_error_EAST.m" and
%                     "get_EFIT_parameters_EAST.m".  These routines were
%                     written by Wang Bo.
%  2016/06/29 - RSG; Add the prad, pinput, and radfrac data.  We had to
%                     reprocess the prad data to avoid excessive and
%                     non-causal smoothing
%  2016/12/19 - RSG; Resample the flattop of all shots at higher rate.
%  2017/04/06 - RSG; Add n=1 data from saddle coils and the RMP coils.
%  2017/12/05 - RSG; Add density and Greenwald fraction
%  2017/12/08 - RSG; Add Wmhd
%  2017/12/11 - RSG; Do not write records to database that have
%                    time_until_disrupt < 0;
%  2017/12/11 - RSG; Added p_ecrh and p_nbi.  Re-ran all shots to populate
%                    these two new parameters.  Also repopulate radiated power
%                    fractions with the baseline offset removed from p_ecrh
%  2018/08/08 - RSG; Added zcur_lmsz and z_error_lmsz.  Re-ran all shots to
%                    populate these two new parameters.
%  2018/10/12 - RSG; Added n_equal_1_normalized and Te_width_normalized
%  2018/12/13 - RSG; Added many new columns for real time (_RT) data from
%                    the EAST PCS.  Also changed the name of the database
%                    we connect to to "east_disruption" after Wang Feng
%                    renamed it.
%  2019/05/10 - RSG; Added new columns for normalized quantities, such as
%                    ip_error_normalized, etc.
%  2019/05/13 - CR;  Add section for Pefit
%  2019/05/14 - RSG; Add section for prad_peaking.  Also, remove the
%                    following columns from the database, and from this
%                    program: dip_dt, dbetap_dt, dli_dt, dWmhd_dt, dn_dt,
%                    and z_times_vz
%  2019/11/06 - RSG; Add section for Mirnov_std and Mirnov_std_normalized,
%                    and load these data for all shots.  Also rerun the
%                    n_equal_1_normalized data for shots > 65165 to account
%                    for Btor "\it" signal changing back to eng_tree.
%  2019/12/23 - RSG; Re-enable the EFIT upper and lower gaps.  Jinxiang Zhu
%                    has written a routine to explicitly calculate these
%                    parameters, since the EAST version of EFIT does not
%                    output them directly.
%  2020/03/11 - CR;  Modified to add kappa_area to EAST SQL db.
%  2021/12/08 - CR;  Populate the database with shots from 2018-2020
%                    Also, add n1rms, n2rms (routines written by J. Zhu)

% Define path to some required Matlab routines

oldpath = addpath('/home/ASIPP/granetz/matlab', ...
                  '/home/ASIPP/granetz/disruption_warning_database', ...
                  '/home/ASIPP/granetz/disruption_database', ...
                  '/home/ASIPP/jinxiang.zhu');
%                 '/home/ASIPP/xianggu/Disruption/code');

% Turn off annoying warning messages about calculations with NaN values in
% the Matlab INTERP1 interpolation routine.

warning_status = warning('query','MATLAB:interp1:NaNinY'); 
if (strcmp(warning_status.state, 'on'));
  warning('off', 'MATLAB:interp1:NaNinY');
end;

% Connect to the SQL database

db = set_database('east_disruption'); % Connect to our EAST SQL database

% Before starting the main loop over the list of shots, get the list of all
% disruption shots and times from the disruption database.  These data will
% be needed to set up non-uniform sampling times prior to disruptions, and
% for filling in the field 'time_until_disrupt' for those shots that
% disrupted.

result = fetch(db, 'select shot,t_disrupt from disruptions order by shot');
disruption_shotlist = cell2mat(result(:,1));
t_disrupt_list = cell2mat(result(:,2));

mdsconnect('mds.ipp.ac.cn'); % Connect to EAST MDSplus server

% A few more minor things to do before we start the main loop over the list
% of shots

nshots = length(shotlist);
ndigits = int2str(floor(log10(nshots)) + 1); % (for formatting in the 
                                             %  upcoming fprintf statement)
%*********************
%
% Start the main loop over the list of shots.
%
%*********************

for ishot = 1:nshots;
  shot = double(shotlist(ishot)); % mdsopen at EAST needs non-integer shot #
  fprintf(1,['Processing shot %7i  (%' ndigits 'i/%' ndigits ...
    'i  %6.2f%%)\n'], shot, ishot, nshots, ishot/nshots*100);

% First, call the routine "check_for_valid_plasma" to see if the plasma
% current exists, and if it satisfies minimum constraints on amplitude and
% duration.  If not, then skip this shot.

  minimum_ip = 200e3;      % Plasma current must be >= 200 kA
  minimum_duration = 0.6;  % Plasma duration must be >= 0.6 s

  [valid, ipmax, duration] = check_for_valid_plasma(shot, minimum_ip, ...
    minimum_duration);

  if (~valid); continue; end;  % If not a valid shot, skip it

% Define the times for the database time slices.  Sampling during the
% rampup, flattop, and rampdown periods can be different, if desired.
% Much faster sampling is done prior to disruptions.

% Call routine "get_flattop_times" to find the start and stop times of the
% flattop.

  [flattop_start, flattop_stop] = get_flattop_times(shot);

% Set the sampling frequency during rampup, flattop, and rampdown.  (These
% may need to be adjusted after studying the data.)

  dt_rampup   = 0.1;  % 0.1 s sampling during rampup
  dt_flattop  = 0.1;  % 0.1 s sampling during flattop
  dt_rampdown = 0.1;  % 0.1 s sampling during rampdown

  if (flattop_start > 0.2);
    times_rampup = [0.2 : dt_rampup : flattop_start];
  else;
    times_rampup = [];
  end;

  if (flattop_stop - flattop_start > 0.0);
    flattop_sampling_start = ceil(flattop_start/dt_flattop) * dt_flattop;
    times_flattop = [flattop_sampling_start : dt_flattop : flattop_stop];
  else;
    times_flattop = [];
  end;

  if (flattop_stop < duration);
    rampdown_sampling_start = ceil(flattop_stop/dt_rampdown) * dt_rampdown;
    times_rampdown = [rampdown_sampling_start : dt_rampdown : duration];
  else;
    times_rampdown = [];
  end;

  times_for_db = unique([times_rampup, times_flattop, times_rampdown]);

% If this shot disrupted, add in additional sampling times prior to the
% disruption at a much higher sampling rate.  (The pre-disruption sampling
% rate and duration may need to be adjusted after studying the data.)

  dt_before_disruption = 0.010;       % Use 10 ms sampling starting
  duration_before_disruption = 0.25;  % 0.25 s before the disruption time

  [disrupted, indx] = ismember(shot, disruption_shotlist);
  if (disrupted);
    tdis = t_disrupt_list(indx);
    times_before_disruption = ...
      [(tdis - duration_before_disruption) : dt_before_disruption : tdis];
    indices = find(times_for_db < min(times_before_disruption));
    times_for_db = unique([times_for_db(indices), times_before_disruption]);
  end;

  ntimes_for_db = length(times_for_db);
  if (ntimes_for_db == 0); continue; end;

% Get list of time slices that are already in the database for this shot
%% CR 2021/12/8 this part needs to change?
  result = fetch(db, ['select time,dbkey from disruption_warning where ' ...
    'shot = ' int2str(shot) ' order by time']);
  if (length(result) > 0);
    times_in_db = cell2mat(result(:,1));
    dbkey = cell2mat(result(:,2));
  else;
    times_in_db = [];
  end;
  ntimes_in_db = length(times_in_db);

% The time slices (times_for_db) need to be separated into two groups:
%   (1) times that already exist in the database (times_in_db), and
%   (2) times that do not exist yet in the database.
% For the 1st group, the existing records will be updated in the database.
% For the 2nd group, new records will be inserted into the database.
% ('Update' and 'insert' are different SQL operations.)  Any time slice
% that is within 10 microseconds of an existing time slice will be
% considered as already existing in the database.

  record_keys = NaN(ntimes_for_db,1);
  if (ntimes_in_db > 0);
    for itime = 1:ntimes_for_db;
      indx = find(abs(times_in_db - times_for_db(itime)) <= 10e-6);
      if ~isempty(indx);
        record_keys(itime) = dbkey(indx);
      end;
    end;
  end;

% Create two cell arrays, 'fields' and 'values'.  The 'fields' array will
% contain the names of the database fields that have data to be written to
% the database.  The 'values' array contains the values for each of the
% fields.

  nfields = length(get_columns(db, 'disruption_warning'));
  fields = cell(1, nfields);
  values = cell(ntimes_for_db, nfields);

% Now begin filling in the two cell arrays, 'fields' and 'values', with all
% the data that are available for this shot.  Robert Granetz and Wang Bo
% wrote most of the Matlab routines that get these data.

  field_counter = 0;

% Start with the time values for the 'time' field.  All of the other data
% will be on this timebase.

    field_counter = field_counter + 1;
    fields(field_counter) = {'time'};
    values(:,field_counter) = num2cell(times_for_db);

% Next, fill in the 'shot' field.

  shot_values = int32(zeros(ntimes_for_db,1)) + shot;

    field_counter = field_counter + 1;
    fields(field_counter) = {'shot'};
    values(:,field_counter) = num2cell(shot_values);

% Next, if this shot disrupted, calculate the values for
% 'time_until_disrupt'.  If the shot did not disrupt, then set the values
% to NaN (Matlab's "Not-a-Number").

  [disrupted, indx] = ismember(shot, disruption_shotlist);
  if (disrupted);
    time_until_disrupt =  t_disrupt_list(indx) - times_for_db;
  else;
    time_until_disrupt = NaN(ntimes_for_db,1);
  end;

    field_counter = field_counter + 1;
    fields(field_counter) = {'time_until_disrupt'};
    values(:,field_counter) = num2cell(time_until_disrupt);

% Next, get the measured plasma current (ip), the programmed plasma current
% (ip_prog), the error and the normalized error between the measured and
% programmed currents (ip_error, ip_error_normalized), and the time
% derivatives of ip and ip_prog (dip_dt and dipprog_dt).  The ip_prog data
% is not included in the database.

  [ip, ip_prog, ip_error, ip_error_normalized, dip_dt, dipprog_dt] = ...
    get_Ip_parameters(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'ip'};
    values(:,field_counter) = num2cell(ip);

    field_counter = field_counter + 1;
    fields(field_counter) = {'ip_error'};
    values(:,field_counter) = num2cell(ip_error);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'ip_error_normalized'};
    values(:,field_counter) = num2cell(ip_error_normalized);

%   field_counter = field_counter + 1;
%   fields(field_counter) = {'dip_dt'};
%   values(:,field_counter) = num2cell(dip_dt);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'dipprog_dt'};
    values(:,field_counter) = num2cell(dipprog_dt);
  
% Get the loop voltage

  v_loop = get_v_loop_EAST(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'v_loop'};
    values(:,field_counter) = num2cell(v_loop);

% Get the Z-coordinate of the plasma centroid (zcur, zcur_lmsz,
% zcur_normalized, zcur_lmsz_normalized), and the error between zcur and
% the programmed zcur (z_error, z_error_lmsz, z_error_normalized,
% z_error_lmsz_normalized)

  [z_error, z_prog, zcur, zcur_lmsz, z_error_lmsz, ...
    zcur_lmsz_normalized, z_error_lmsz_normalized] = ...
    get_Z_error_EAST(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'zcur'};
    values(:,field_counter) = num2cell(zcur);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'z_error'};
    values(:,field_counter) = num2cell(z_error);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'zcur_lmsz'};
    values(:,field_counter) = num2cell(zcur_lmsz);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'z_error_lmsz'};
    values(:,field_counter) = num2cell(z_error_lmsz);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'zcur_lmsz_normalized'};
    values(:,field_counter) = num2cell(zcur_lmsz_normalized);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'z_error_lmsz_normalized'};
    values(:,field_counter) = num2cell(z_error_lmsz_normalized);
  
% Get vertical drift velocity and associated parameter

%  [Z_cur1, v_z, z_times_v_z] = get_Z(shot, times_for_db);

%    field_counter = field_counter + 1;
%    fields(field_counter) = {'v_z'};
%    values(:,field_counter) = num2cell(v_z);
  
%   field_counter = field_counter + 1;
%   fields(field_counter) = {'z_times_v_z'};
%   values(:,field_counter) = num2cell(z_times_v_z);

% Get density and Greenwald fraction

  [ne, Greenwald_fraction, dn_dt] = ...
    get_density_parameters_EAST(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'n_e'};
    values(:,field_counter) = num2cell(ne);

    field_counter = field_counter + 1;
    fields(field_counter) = {'Greenwald_fraction'};
    values(:,field_counter) = num2cell(Greenwald_fraction);

%   field_counter = field_counter + 1;
%   fields(field_counter) = {'dn_dt'};
%   values(:,field_counter) = num2cell(dn_dt);

% Get heating input powers and radiated output power

  [p_rad, p_ecrh, p_lh, p_oh, p_icrf, p_nbi, rad_input_frac, ...
    rad_loss_frac, p_input] = get_power_EAST(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'p_rad'};
    values(:,field_counter) = num2cell(p_rad);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'p_lh'};
    values(:,field_counter) = num2cell(p_lh);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'p_oh'};
    values(:,field_counter) = num2cell(p_oh);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'p_icrf'};
    values(:,field_counter) = num2cell(p_icrf);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'p_ecrh'};
    values(:,field_counter) = num2cell(p_ecrh);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'p_nbi'};
    values(:,field_counter) = num2cell(p_nbi);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'rad_input_frac'};
    values(:,field_counter) = num2cell(rad_input_frac);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'rad_loss_frac'};
    values(:,field_counter) = num2cell(rad_loss_frac);

% Get n=1 error field data

  [n_equal_1_normalized, n_equal_1_mode, n_equal_1_phase, rmp_n_equal_1, ...
    rmp_n_equal_1_phase, btor] = get_n_equal_1_data_EAST(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'n_equal_1_normalized'};
    values(:,field_counter) = num2cell(n_equal_1_normalized);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'n_equal_1_mode'};
    values(:,field_counter) = num2cell(n_equal_1_mode);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'n_equal_1_phase'};
    values(:,field_counter) = num2cell(n_equal_1_phase);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'rmp_n_equal_1'};
    values(:,field_counter) = num2cell(rmp_n_equal_1);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'rmp_n_equal_1_phase'};
    values(:,field_counter) = num2cell(rmp_n_equal_1_phase);

    field_counter = field_counter + 1;
    fields(field_counter) = {'btor'};
    values(:,field_counter) = num2cell(btor);

% Get PEFIT parameters

  [beta_n, beta_p, kappa, li, q95, Wmhd, aminor] = ...
    get_PEFIT_parameters(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'Pbeta_n'};
    values(:,field_counter) = num2cell(beta_n);
    
    field_counter = field_counter + 1;
    fields(field_counter) = {'Pbeta_p'};
    values(:,field_counter) = num2cell(beta_p);
    
    field_counter = field_counter + 1;
    fields(field_counter) = {'Pkappa'};
    values(:,field_counter) = num2cell(kappa);
    
    field_counter = field_counter + 1;
    fields(field_counter) = {'Pli'};
    values(:,field_counter) = num2cell(li);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'Pq95'};
    values(:,field_counter) = num2cell(q95);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'PWmhd'};
    values(:,field_counter) = num2cell(Wmhd);

    field_counter = field_counter + 1;
    fields(field_counter) = {'Paminor'};
    values(:,field_counter) = num2cell(aminor);


% Get EFIT parameters

  [beta_n, beta_p, dbetap_dt, kappa, upper_gap, lower_gap, ...
    li, dli_dt, q0, qstar, q95, Wmhd, dWmhd_dt] = ...
    get_EFIT_parameters_EAST(shot, times_for_db);


    field_counter = field_counter + 1;
    fields(field_counter) = {'beta_n'};
    values(:,field_counter) = num2cell(beta_n);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'beta_p'};
    values(:,field_counter) = num2cell(beta_p);
  
%   field_counter = field_counter + 1;
%   fields(field_counter) = {'dbetap_dt'};
%   values(:,field_counter) = num2cell(dbetap_dt);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'kappa'};
    values(:,field_counter) = num2cell(kappa);
  


    field_counter = field_counter + 1;
    fields(field_counter) = {'upper_gap'};
    values(:,field_counter) = num2cell(upper_gap);

    field_counter = field_counter + 1;
    fields(field_counter) = {'lower_gap'};
    values(:,field_counter) = num2cell(lower_gap);
  


    field_counter = field_counter + 1;
    fields(field_counter) = {'li'};
    values(:,field_counter) = num2cell(li);
  
%   field_counter = field_counter + 1;
%   fields(field_counter) = {'dli_dt'};
%   values(:,field_counter) = num2cell(dli_dt);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'q0'};
    values(:,field_counter) = num2cell(q0);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'qstar'};
    values(:,field_counter) = num2cell(qstar);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'q95'};
    values(:,field_counter) = num2cell(q95);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'Wmhd'};
    values(:,field_counter) = num2cell(Wmhd);

%   field_counter = field_counter + 1;
%   fields(field_counter) = {'dWmhd_dt'};
%   values(:,field_counter) = num2cell(dWmhd_dt);


% Get kappa_area
  [kappa_area] = get_kappa_area(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'kappa_area'};
    values(:,field_counter) = num2cell(kappa_area);

%{
% Get toroidal rotation velocities (on-axis and mid-radius)

  [v_0, v_mid] = get_rotation_velocity(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'v_0'};
    values(:,field_counter) = num2cell(v_0);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'v_mid'};
    values(:,field_counter) = num2cell(v_mid);
  
%}
  [Te_width_normalized, Te_width_MI] = get_Te_width(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'Te_width_normalized'};
    values(:,field_counter) = num2cell(Te_width_normalized);

    field_counter = field_counter + 1;
    fields(field_counter) = {'Te_width_MI'};
    values(:,field_counter) = num2cell(Te_width_MI);

  [p_rad_RT, p_ecrh_RT, p_lh_RT, p_oh_RT, p_icrf_RT, p_nbi_RT, ...
    rad_input_frac_RT, rad_loss_frac_RT, ip_error_RT, ...
    q95_RT, beta_p_RT, li_RT, Wmhd_RT] = get_PCS_EAST(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'p_rad_RT'};
    values(:,field_counter) = num2cell(p_rad_RT);

    field_counter = field_counter + 1;
    fields(field_counter) = {'ip_error_RT'};
    values(:,field_counter) = num2cell(ip_error_RT);

    field_counter = field_counter + 1;
    fields(field_counter) = {'q95_RT'};
    values(:,field_counter) = num2cell(q95_RT);

    field_counter = field_counter + 1;
    fields(field_counter) = {'beta_p_RT'};
    values(:,field_counter) = num2cell(beta_p_RT);

    field_counter = field_counter + 1;
    fields(field_counter) = {'li_RT'};
    values(:,field_counter) = num2cell(li_RT);

    field_counter = field_counter + 1;
    fields(field_counter) = {'Wmhd_RT'};
    values(:,field_counter) = num2cell(Wmhd_RT);

    field_counter = field_counter + 1;
    fields(field_counter) = {'p_nbi_RT'};
    values(:,field_counter) = num2cell(p_nbi_RT);

    field_counter = field_counter + 1;
    fields(field_counter) = {'p_lh_RT'};
    values(:,field_counter) = num2cell(p_lh_RT);

  prad_peaking = get_prad_peaking(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'prad_peaking'};
    values(:,field_counter) = num2cell(prad_peaking);

  [Mirnov_std_normalized, Mirnov_std] = ...
    get_Mirnov_std(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'Mirnov_std'};
    values(:,field_counter) = num2cell(Mirnov_std);

    field_counter = field_counter + 1;
    fields(field_counter) = {'Mirnov_std_normalized'};
    values(:,field_counter) = num2cell(Mirnov_std_normalized);

  [n1rms, n2rms, n1rms_normalized, n2rms_normalized] = ...
    get_n1rms_n2rms_EAST(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'n1rms'};
    values(:,field_counter) = num2cell(Mirnov_std);

    field_counter = field_counter + 1;
    fields(field_counter) = {'n1rms_normalized'};
    values(:,field_counter) = num2cell(Mirnov_std_normalized);

    field_counter = field_counter + 1;
    fields(field_counter) = {'n2rms'};
    values(:,field_counter) = num2cell(Mirnov_std);

    field_counter = field_counter + 1;
    fields(field_counter) = {'n2rms_normalized'};
    values(:,field_counter) = num2cell(Mirnov_std_normalized);

  H98 = get_H98_EAST(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'H98'};
    values(:,field_counter) = num2cell(prad_peaking);

% All the information for this shot has been obtained.  Now prepare to
% enter the information into the disruption_warning table.

  if (field_counter == 0); % If no data is available for any of the fields,
    continue;              % then skip to the next shot in the list
  end;

  fields = fields(1:field_counter);   % Truncate these cell arrays to get
  values = values(:,1:field_counter); % rid of unassigned fields.

% But first, deal with rare problem when some values can be infinite, or
% can exceed the maximum allowed value for MySQL float type (~1e38)

  for ifield = 1:length(fields);
    dummy = cell2mat(values(:,ifield));
    bad_vals = find(abs(dummy) > 1.e30);
    dummy(bad_vals) = 1.e30 .* sign(dummy(bad_vals));
    values(:,ifield) = num2cell(dummy);
  end;

% For the time slices that already exist in the database, an SQL 'update'
% operation will be performed.  Otherwise, an SQL 'insert' operation will
% be performed.

  for itime = 1:ntimes_for_db;
    if (isnan(record_keys(itime)));
      insert(db, 'disruption_warning', fields, values(itime,:));
    else;

% Handling of NaN by the Matlab database 'update' routine changed in 2016.
% The NaN values have to be changed to an empty value ''.

%     null_indices = find(isnan(cell2mat(values(itime,:))));
%     values(itime, null_indices) = {''};

      for iparam = 1:length(values(itime,:));
        if isnan(cell2mat(values(itime, iparam)));
          values(itime, iparam) = {''};
        end;
      end;

      update(db, 'disruption_warning', fields, values(itime,:), ...
        ['where dbkey = ' num2str(record_keys(itime),'%i')]);
    end;
  end;

% Okay, this shot is done.  Go to the next shot in shotlist.

end;

% Okay, we have processed all the shots in the list.  We just need to clean
% up a few things, and then we're done.

path(oldpath);
close(db);
mdsdisconnect;
if (strcmp(warning_status.state, 'on'));  % Turn warning messages about calcs
  warning('on', 'MATLAB:interp1:NaNinY'); % w/NaN's back on if it was on
end;                                      % when this routine was called.
