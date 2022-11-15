function disruption_warning_database_d3d_montes(shotlist);

% 'shotlist' is a Matlab array of shot numbers.  It can be a single shot,
%   or a sequence of shots [160946:161837], or a list of specific shots
%   [160988, 161129, 165492], or combinations of these [160946:161837, 160988,
%   161129, 165492].
%
% Examples of usage:
%
%  disruption_warning_database_d3d(160946);  % Process shot 160946 only
%
%  disruption_warning_database_d3d(160946:161837);  % Process all plasma
%                                                   % shots sequentially
%                                                   % from 160946 to 161837
%
%  disruption_warning_database_d3d(shotlist); % Process the list of shots in
%                                             % the Matlab array 'shotlist'
%
% Revision/history:
%  2016/03    - RSG; Created main routine (adapted from the disruption
%                     warning database programs that I wrote for the C-Mod
%                     and EAST tokamaks)
%  2017/05/02 - RSG; Added Te_peaking parameter

% Define path on the Iris cluster to some required Matlab routines

%oldpath = addpath( ...
%  '/home/granetzr/matlab', ...
oldpath = addpath('/fusion/projects/disruption_warning/matlab_programs');

% Turn off annoying warning messages about calculations with NaN values in
% the Matlab INTERP1 interpolation routine.

%{
warning_status = warning('query','MATLAB:interp1:NaNinY'); 
if (strcmp(warning_status.state, 'on'));
  warning('off', 'MATLAB:interp1:NaNinY');
end;
%}

warning_status = warning;
warning('off')

% Before starting the main loop over the list of shots, get the list of all
% disruption shots and times from the disruption database.  These data will
% be needed to set up more frequent sampling times prior to disruptions
% and for filling in the field 'time_until_disrupt' for those shots that
% disrupted.

db = set_database('d3drdb'); % Connect to the SQL database that contains
                             % the 'disruptions' table
result = fetch(db, 'select shot, t_disrupt from disruptions order by shot');
close(db);
disruption_shotlist = int32(cell2mat(result(:,1)));
t_disrupt_list = cell2mat(result(:,2));

% Connect to the SQL database that contains the 'disruption_warning' table
%-----------
% For testing purposes only, we can connect instead to the 'test_db' SQL
% database and use the 'test_disruption_warning' table
%-----------

db = set_database('d3drdb');
%db = set_database('test_db');

mdsconnect('atlas.gat.com'); % Connect to DIII-D MDSplus server

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
  shot = shotlist(ishot);
  fprintf(1,['Processing shot %7i  (%' ndigits 'i/%' ndigits ...
    'i  %6.2f%%)\n'], shot, ishot, nshots, ishot/nshots*100);

% First, call the routine "check_for_valid_plasma" to see if the plasma
% current exists, and if it satisfies minimum constraints on amplitude and
% duration.  If not, then skip this shot.
% Specifications for minimum pulse_length and minimum ip were decided in a
% conversation I had with Nick Eidietis on 2016/02/25.

  minimum_ip = 400e3;      % Plasma current must be >= 400 kA
  minimum_duration = 0.1;  % Plasma duration must be >= 0.1 s

  [validity, ipmax, duration] = check_for_valid_plasma(shot, minimum_ip, ...
    minimum_duration);

  if (~validity); continue; end;  % If not a valid shot, skip it

% Define the times for the database time slices.  Create a timebase with a
% suitable sampling rate, and then add in more frequent time slices prior
% to disruptions.  After studying the initial dataset, we may decide to
% change the sampling parameters, and repopulate the database.

  times_for_db = [0.100 : 0.025 : duration]; % Use 25 ms sampling for now

% If this shot disrupted, add in additional sampling times prior to the
% disruption at a much higher sampling rate.  (The pre-disruption sampling
% rate and duration may need to be adjusted after studying the data.)

  dt_before_disruption = 0.002;       % Use 2 ms sampling starting at
  duration_before_disruption = 0.10;  % 0.10 s before the disruption time

  [disrupted, indx] = ismember(shot, disruption_shotlist);
  if (disrupted);
    tdis = t_disrupt_list(indx);
    times_before_disruption = ...
      [(tdis - duration_before_disruption) : dt_before_disruption : tdis];
%  We have decided that we do not want any overlap between the standard
%  sample times and the pre-disruption sample times.  Therefore, we remove
%  all standard sample times that occur during the pre-disruption period.
    okay_indices = find(times_for_db < (tdis - duration_before_disruption));
    times_for_db = [times_for_db(okay_indices), times_before_disruption];
  end;

  ntimes_for_db = length(times_for_db);
  if (ntimes_for_db == 0); continue; end;

% Get list of time slices that are already in the database for this shot

  result = fetch(db, ['select time, dbkey from disruption_warning ' ...
    'where shot = ' int2str(shot) ' order by time']);
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
% ('update' and 'insert' are different SQL operations.)  Any time slice
% that is within 10 microseconds of an existing time slice will be
% considered as already existing in the database.

  record_keys = NaN(ntimes_for_db, 1);
  if (ntimes_in_db > 0);
    for itime = 1:ntimes_for_db;
      indx = find(abs(times_in_db - times_for_db(itime)) <= 10e-6);
      if ~isempty(indx);
        record_keys(itime) = dbkey(indx(1));
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
% the data that are available for this shot.  Robert Granetz, Alex
% Tinguely, and Cristina Rea wrote most of the Matlab routines that get
% these data.

  field_counter = 0;
%{
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
%}
% Next, get the measured plasma current (ip), the programmed plasma current
% (ip_prog), the error between the measured and programmed currents
% (ip_error), and the time derivatives of ip and ip_prog (dip_dt and
% dipprog_dt).  The ip_prog data is not included in the database.
%{
  [ip, ip_prog, ip_error, dip_dt, dipprog_dt, power_supply_railed] = ...
    get_Ip_parameters(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'ip'};
    values(:,field_counter) = num2cell(ip);

    field_counter = field_counter + 1;
    fields(field_counter) = {'ip_error'};
    values(:,field_counter) = num2cell(ip_error);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'dip_dt'};
    values(:,field_counter) = num2cell(dip_dt);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'dipprog_dt'};
    values(:,field_counter) = num2cell(dipprog_dt);

    field_counter = field_counter + 1;
    fields(field_counter) = {'power_supply_railed'};
    values(:,field_counter) = num2cell(power_supply_railed);

% Get the Z-coordinate of the plasma centroid (zcur), and the error between
% zcur and the programmed zcur (z_error)

  [z_error, z_prog, zcur] = get_Z_error_d3d(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'zcur'};
    values(:,field_counter) = num2cell(zcur);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'z_error'};
    values(:,field_counter) = num2cell(z_error);

% Get n=1 locked mode data

% [n_equal_1_mode, IRLM] = get_locked_mode_data_d3d(shot, times_for_db);

  [n_equal_1_mode, n_equal_1_normalized] = ...
    get_n1_bradial_d3d(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'n_equal_1_mode'};
    values(:,field_counter) = num2cell(n_equal_1_mode);

    field_counter = field_counter + 1;
    fields(field_counter) = {'n_equal_1_normalized'};
    values(:,field_counter) = num2cell(n_equal_1_normalized);

% Get vertical drift velocity and associated parameter

  [Z_cur1, v_z, z_times_v_z] = get_Z(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'v_z'};
    values(:,field_counter) = num2cell(v_z);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'z_times_v_z'};
        field_counter = field_counter + 1;
    fields(field_counter) = {'n_equal_1_mode'};
    values(:,field_counter) = num2cell(n_equal_1_mode);


values(:,field_counter) = num2cell(z_times_v_z);

% Get heating input powers and radiated output power

  [p_rad, p_nbi, p_ohm, p_ech, radiated_fraction, p_input, v_loop] = ...
    get_power_d3d(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'p_rad'};
    values(:,field_counter) = num2cell(p_rad);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'p_nbi'};
    values(:,field_counter) = num2cell(p_nbi);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'p_ohm'};
    values(:,field_counter) = num2cell(p_ohm);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'p_ech'};
    values(:,field_counter) = num2cell(p_ech);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'radiated_fraction'};
    values(:,field_counter) = num2cell(radiated_fraction);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'v_loop'};
    values(:,field_counter) = num2cell(v_loop);
 
% Get EFIT parameters

  [beta_n, beta_p, dbetap_dt, kappa, upper_gap, lower_gap, ...
    li, dli_dt, q0, qstar, q95, v_loop_efit, Wmhd, dWmhd_dt] = ...
    get_EFIT_parameters(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'beta_n'};
    values(:,field_counter) = num2cell(beta_n);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'beta_p'};
    values(:,field_counter) = num2cell(beta_p);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'dbetap_dt'};
    values(:,field_counter) = num2cell(dbetap_dt);
  
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
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'dli_dt'};
    values(:,field_counter) = num2cell(dli_dt);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'q0'};
    values(:,field_counter) = num2cell(q0);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'qstar'};
    values(:,field_counter) = num2cell(qstar);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'q95'};
    values(:,field_counter) = num2cell(q95);
  
%   field_counter = field_counter + 1;
%   fields(field_counter) = {'v_loop_efit'};
%   values(:,field_counter) = num2cell(v_loop_efit);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'Wmhd'};
    values(:,field_counter) = num2cell(Wmhd);

    field_counter = field_counter + 1;
    fields(field_counter) = {'dWmhd_dt'};
    values(:,field_counter) = num2cell(dWmhd_dt);

% Get toroidal rotation velocities (on-axis and mid-radius)

  [v_0, v_mid] = get_rotation_velocity(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'v_0'};
    values(:,field_counter) = num2cell(v_0);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'v_mid'};
    values(:,field_counter) = num2cell(v_mid);
  
% Get density and Greenwald fraction

  [ne, Greenwald_fraction, dndt] = get_density_parameters(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'n_e'};
    values(:,field_counter) = num2cell(ne);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'Greenwald_fraction'};
    values(:,field_counter) = num2cell(Greenwald_fraction);

    field_counter = field_counter + 1;
    fields(field_counter) = {'dn_dt'};
    values(:,field_counter) = num2cell(dndt);

% Get Te profile width from 'blessed' Thomson scattering data in the MDS
% electrons tree

  [Te_HWHM, Te_width_normalized] = get_TS_data(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'Te_HWHM'};
    values(:,field_counter) = num2cell(Te_HWHM);

    field_counter = field_counter + 1;
    fields(field_counter) = {'Te_width_normalized'};
    values(:,field_counter) = num2cell(Te_width_normalized);
%}
% Get the measured plasma current, the programmed plasma current, the error
% between the measured and programmed currents, and the time derivatives of
% ip and ip_prog from real-time data available in the PCS.  The PCS data
% are what we will actually use in our real-time algorithm, and this could
% be slightly different than the data in MDS.  The
%{
  [ip_RT, ip_prog_RT, ip_error_RT, dip_dt_RT, dipprog_dt_RT, ...
    power_supply_railed] = get_Ip_parameters_RT(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'ip_RT'};
    values(:,field_counter) = num2cell(ip_RT);

    field_counter = field_counter + 1;
    fields(field_counter) = {'ip_error_RT'};
    values(:,field_counter) = num2cell(ip_error_RT);
  
%   field_counter = field_counter + 1;
%   fields(field_counter) = {'dip_dt_RT'};
%   values(:,field_counter) = num2cell(dip_dt);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'dipprog_dt_RT'};
    values(:,field_counter) = num2cell(dipprog_dt_RT);

%   field_counter = field_counter + 1;
%   fields(field_counter) = {'power_supply_railed'};
%   values(:,field_counter) = num2cell(power_supply_railed);

% Get density and Greenwald fraction from real-time data available in the
% PCS.  The PCS data are what we will actually use in our real-time
% algorithm, and this could be slightly different than the data in MDS.

  [ne_RT, Greenwald_fraction_RT, dndt_RT] = ...
    get_density_parameters_RT(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'n_e_RT'};
    values(:,field_counter) = num2cell(ne_RT);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'Greenwald_fraction_RT'};
    values(:,field_counter) = num2cell(Greenwald_fraction_RT);

%   field_counter = field_counter + 1;
%   fields(field_counter) = {'dn_dt_RT'};
%   values(:,field_counter) = num2cell(dndt_RT);

% Get Te profile width from real-time Thomson data available in the PCS.
% The PCS data are what we will actually use in our real-time algorithm,
% and this could be slightly different than the 'blessed' data in MDS.
% (But we have compared these on a subset of shots and not found
% significant differences.)

  [Te_HWHM_RT, Te_width_normalized_RT] = get_TS_data_RT(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'Te_HWHM_RT'};
    values(:,field_counter) = num2cell(Te_HWHM_RT);

    field_counter = field_counter + 1;
    fields(field_counter) = {'Te_width_normalized_RT'};
    values(:,field_counter) = num2cell(Te_width_normalized_RT);
    
% Get EFIT_RT parameters

  [beta_n_RT, beta_p_RT, dbetap_dt_RT, kappa_RT, upper_gap_RT, ...
    lower_gap_RT, li_RT, dli_dt_RT, q0_RT, qstar_RT, q95_RT, ...
    v_loop_efit_RT, Wmhd_RT, dWmhd_dt_RT] = ...
    get_EFIT_parameters_RT(shot, times_for_db);

%   field_counter = field_counter + 1;
%   fields(field_counter) = {'beta_n_RT'};
%   values(:,field_counter) = num2cell(beta_n_RT);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'beta_p_RT'};
    values(:,field_counter) = num2cell(beta_p_RT);
  
%   field_counter = field_counter + 1;
%   fields(field_counter) = {'dbetap_dt_RT'};
%   values(:,field_counter) = num2cell(dbetap_dt_RT);
  
%   field_counter = field_counter + 1;
%   fields(field_counter) = {'kappa_RT'};
%   values(:,field_counter) = num2cell(kappa_RT);
  
%   field_counter = field_counter + 1;
%   fields(field_counter) = {'upper_gap_RT'};
%   values(:,field_counter) = num2cell(upper_gap_RT);

%   field_counter = field_counter + 1;
%   fields(field_counter) = {'lower_gap_RT'};
%   values(:,field_counter) = num2cell(lower_gap_RT);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'li_RT'};
    values(:,field_counter) = num2cell(li_RT);
  
%   field_counter = field_counter + 1;
%   fields(field_counter) = {'dli_dt_RT'};
%   values(:,field_counter) = num2cell(dli_dt_RT);
  
%   field_counter = field_counter + 1;
%   fields(field_counter) = {'q0_RT'};
%   values(:,field_counter) = num2cell(q0_RT);
  
%   field_counter = field_counter + 1;
%   fields(field_counter) = {'qstar_RT'};
%   values(:,field_counter) = num2cell(qstar_RT);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'q95_RT'};
    values(:,field_counter) = num2cell(q95_RT);
  
%   field_counter = field_counter + 1;
%   fields(field_counter) = {'v_loop_efit_RT'};
%   values(:,field_counter) = num2cell(v_loop_efit_RT);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'Wmhd_RT'};
    values(:,field_counter) = num2cell(Wmhd_RT);

%   field_counter = field_counter + 1;
%   fields(field_counter) = {'dWmhd_dt_RT'};
%   values(:,field_counter) = num2cell(dWmhd_dt_RT);

%}

% Grab peaking factors for electron density and temperature
  [Te_PF_CVA, ne_PF_CVA, ~] = ...
    get_Thomson_peaking_factors(shot, times_for_db,'core_vs_all');

    field_counter = field_counter + 1;
    fields(field_counter) = {'Te_peaking_CVA'};
    values(:,field_counter) = num2cell(Te_PF_CVA);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'ne_peaking_CVA'};
    values(:,field_counter) = num2cell(ne_PF_CVA);
%{
    field_counter = field_counter + 1;
    fields(field_counter) = {'Te_peaking_CVE'};
    values(:,field_counter) = num2cell(Te_PF_CVE);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'ne_peaking_CVE'};
    values(:,field_counter) = num2cell(ne_PF_CVE);
%}
% All the information for this shot has been obtained.  Now prepare to
% enter the information into the disruption_warning table.

  if (field_counter == 0); % If no data is available for any of the fields,
    continue;              % then skip to the next shot in the list
  end;

  fields = fields(1:field_counter);   % Truncate these cell arrays to get
  values = values(:,1:field_counter); % rid of unassigned fields.

% For the time slices that already exist in the database, an SQL 'update'
% operation will be performed.  Otherwise, an SQL 'insert' operation will
% be performed.

  for itime = 1:ntimes_for_db;
    if (isnan(record_keys(itime)));
      insert(db, 'disruption_warning', fields, values(itime,:));
    else;

% Apparently in the latest version of Matlab's Database Toolbox, the
% "update" function no longer handles 'NaN' values.  These now have to be
% replaced with empty values ('').  Here is the exact statement in the
% "update" function documentation (for R2016a):
%
%   The database preference settings NullNumberWrite and NullStringWrite do
%   not apply to this function. If data contains null entries and NaNs,
%   convert these entries to an empty value ''.
%
% So I'll temporarily insert the following lines to make my code compatible
% with this reduced capability!
%-------------
      for ifield = 1:field_counter;
        if (isnan(cell2mat(values(itime, ifield))) || ...
            isinf(cell2mat(values(itime, ifield))));
          values(itime, ifield) = {''};
        end;
      end;
%-------------
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

%{
if (strcmp(warning_status.state, 'on'));  % Turn warning messages about calcs
  warning('on', 'MATLAB:interp1:NaNinY'); % w/NaN's back on if it was on
end;                                      % when this routine was called.
%}

warning(warning_status);
end
