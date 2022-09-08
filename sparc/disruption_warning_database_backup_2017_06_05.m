function disruption_warning_database(shotlist, varargin);

% 'shotlist' is a Matlab array of shot numbers.  It can be a single shot,
%   or a sequence of shots [1150901001:1150901099], or a list of specific
%   shots [1150922001, 1150922020, 1150923005], or combinations of these
%   [1150901001:1150901099, 1150922001, 1150922020, 1150923005 ].  Shot #0
%   refers to the current shot.
%
% 'varargin' allows for explicit specification of the fields that will be
%   written into the disruption_warning table.  The varargin parameter(s)
%   must be character strings.  If varargin is not specified, then all
%   fields with available data will be written.
%
% Examples of usage:
%
%  disruption_warning_database(0);
%                         % Process current shot.  First check to see if
%                         % there was a plasma.  (There was a plasma if
%                         % EFIT data exists.)  If there was no plasma,
%                         % then skip the shot.  If no records already
%                         % exist in the database with the specified shot
%                         % and times (within 10 microseconds), then
%                         % insert new records.  If records for the
%                         % specified shot and time(s) already exist, then
%                         % update the existing records.  If varargin
%                         % parameter(s) is/are specified, then only write
%                         % values into the specified fields.
%
%  disruption_warning_database(1150901001:1150901099);
%                         % Process all plasma shots from 1150901001 to
%                         % 1150901099.  
%
%  disruption_warning_database(1150901001:1150901099, 'ip_error',
%    'beta_n', 'z_error');
%                         % Process all plasma shots from 1150901001 to
%                         % 1150901099.  Only write data into fields
%                         % 'ip_error', 'beta_n', and 'z_error'.
%
% Revision/history:
%  2015/09    - RSG; Created (starting from the disruption database program
%                     that I wrote for the EAST tokamak), and using the
%                     many routines that Alex Tinguely wrote to get the
%                     data.
%  2016/01    - RSG; Added 1 ms sampling during the 20 ms period before
%                     disruptions.  This was accomplished by rerunning EFIT
%                     on every disruption, using both the standard timebase
%                     AND the higher sampling timebase just prior to
%                     disruptions.  These EFIT data are written in my EFIT18
%                     tree.  The EFIT18 trees need to be used on disruption
%                     shots, and the ANALYSIS tres need to be used on
%                     non-disruptive shots.
%  2016/02    - RSG; Obtain loop voltage from the mflux calculation, instead
%                     of from EFIT.  This is because EFIT's loop voltage
%                     signal uses a wide acausal filter window which is
%                     unacceptable before disruptions.  This change to the
%                     loop voltage is also incorporated into the get_P_ohm
%                     routine.
%  2016/03/14 - RSG; Added EFIT ssep variable, and changed Wth to Wmhd, and
%                     dWdiam_dt to dWmhd_dt

fieldlist = lower(strtrim(varargin)); % lowercase; remove spaces

% Define path to some required Matlab routines

oldpath = path;
path('/home/granetz/matlab', path);
addpath(genpath('/home/tinguely/Disruptions/Code/Disruption_Database'));
addpath('/home/montes/Disruption_warning_code');

% Turn off annoying warning messages from the "interp1" routine about
% calculations with NaN's.  The original setting for this feature will be
% restored at the end of this program.

warning_status = warning('query','MATLAB:interp1:NaNinY'); 
if (strcmp(warning_status.state, 'on'));
  warning('off', 'MATLAB:interp1:NaNinY');
end;

% Connect to the SQL database and determine how many fields are in the
% disruption_warning table.

db = set_database('logbook'); % Connect to the C-Mod SQL database

fields_in_db = get_columns(db, 'disruption_warning');
fields_in_db = lower(fields_in_db); % lowercase

if (length(fieldlist) == 0);
  nfields = length(fields_in_db);  
else;
  not_okay = ~ismember(fieldlist, fields_in_db);
  if (sum(not_okay) > 0);
    fprintf(1, '\nThese invalid field names will be ignored:\n');
    fprintf(1, '   %s\n', fieldlist{1, not_okay});
    fprintf(1, '\n');
    fieldlist = fieldlist(~not_okay);
  end;
  nfields = length(fieldlist);
end;

% Before starting the main loop over the list of shots, get the list of all
% disruption shots and times from the disruption database.  These data will
% be needed for the field 'time_until_disrupt' for those shots that
% disrupted.

result = fetch(db, 'select shot,t_disrupt from disruptions order by shot');
disruption_shotlist = int32(cell2mat(result(:,1)));
t_disrupt_list = cell2mat(result(:,2));

mdsconnect('alcdata.psfc.mit.edu'); % Connect to C-Mod MDSplus server

% Just a few more things to do before we start the main loop over the list
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
  shot = int32(shotlist(ishot));
  if (shot == 0);
    shot = int32(mdsvalue('current_shot("cmod")'));
  end;

  fprintf(1,['Processing shot %10i  (%' ndigits 'i/%' ndigits ...
    'i  %6.2f%%)\n'], shot, ishot, nshots, ishot/nshots*100);

% First, see if there are any EFIT data for this shot.  If not, then assume
% that there was no significant plasma on this shot, and skip over it.  For
% non-disruptive shots, the EFIT data in the ANALYSIS tree should be used.
% For shots that disrupted, the EFIT data in the EFIT18 tree should be
% used, since it has a burst of fast sampling prior to the disruption time,
% in addition to the standard sampling times throughout the shot.

  if (ismember(shot, disruption_shotlist));
    [shotopened, status] = mdsopen('efit18', shot);
    if (mod(status,2) == 0);
      fprintf(1,'  WARNING: no EFIT18 tree for this disruption shot\n');
      [shotopened, status] = mdsopen('analysis', shot);
    end;
  else;
    [shotopened, status] = mdsopen('analysis', shot);
  end;

  if (mod(status,2) == 1);
    [efittimes, status] = mdsvalue('\efit_aeqdsk:time');
    if (mod(status,2) == 1);
      if (length(efittimes) == 0);
        mdsclose;
        continue;
      end;
    else;
      mdsclose;
      continue;
    end;
  else;
    continue; % If no EFIT tree, skip this shot.
  end;

% Define the times for the database time slices to be the same as the EFIT
% times.  The standard automatic EFIT (in the C-Mod ANALYSIS tree) runs at
% 20 ms intervals, starting at t=0.060 s.  In addition, for disruption
% shots, I have run EFITs (in the EFIT18 tree) at 1 ms intervals, starting
% 20 ms before t_disrupt and stopping at t_disrupt (inclusive), in addition
% to the standard EFIT times.

  times_for_db = efittimes;
  ntimes_for_db = length(times_for_db);

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

  fields = cell(1, nfields);
  values = cell(ntimes_for_db, nfields);

% Now begin filling in the two cell arrays, 'fields' and 'values', with all
% the data that are available for this shot.  Alex Tinguely wrote most of
% the Matlab routines that get these data.  I have modified some of them.

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
% to NaN (Matlab's "Not-a-Number").  Microsoft SqlServer does not support
% the Matlab (and IEEE) infinity value, "Inf".

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
% (ip_prog), the error between the measured and programmed currents
% (ip_error), and the time derivatives of ip and ip_prog (dip_dt and
% dipprog_dt).  The ip_prog data is not included in the database.

  [ip, ip_prog, ip_error, dip_dt, dipprog_dt] = ...
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
  
% Get the Z-coordinate of the plasma centroid (zcur), the programmed zcur
% (Z_prog), the error between zcur and the programmed zcur (z_error), the
% vertical velocity (v_z = d(zcur)/dt), and v_z * zcur

[z_error, Z_prog, zcur, v_z, z_times_v_z] = ...
    get_Z_parameters(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'zcur'};
    values(:,field_counter) = num2cell(zcur);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'z_error'};
    values(:,field_counter) = num2cell(z_error);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'v_z'};
    values(:,field_counter) = num2cell(v_z);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'z_times_v_z'};
    values(:,field_counter) = num2cell(z_times_v_z);
  
% Get heating input powers and radiated output power

  [p_rad, dprad_dt, p_lh, p_oh, p_icrf, p_input, radiated_fraction, ...
    v_loop] = get_power(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'p_rad'};
    values(:,field_counter) = num2cell(p_rad);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'dprad_dt'};
    values(:,field_counter) = num2cell(dprad_dt);
  
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
    fields(field_counter) = {'radiated_fraction'};
    values(:,field_counter) = num2cell(radiated_fraction);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'v_loop'};
    values(:,field_counter) = num2cell(v_loop);

% Get EFIT parameters

  [beta_n, beta_p, dbetap_dt, kappa, upper_gap, lower_gap, ...
    li, dli_dt, q0, qstar, q95, v_loop_efit, Wmhd, dWmhd_dt, ssep] = ...
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
    fields(field_counter) = {'dWmhd_dt'};
    values(:,field_counter) = num2cell(dWmhd_dt);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'ssep'};
    values(:,field_counter) = num2cell(ssep);
  
% Get toroidal rotation velocities (on-axis and mid-radius)

  [v_0, v_mid] = get_rotation_velocity(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'v_0'};
    values(:,field_counter) = num2cell(v_0);

    field_counter = field_counter + 1;
    fields(field_counter) = {'v_0_uncalibrated'};
    values(:,field_counter) = num2cell(v_0);
  
    field_counter = field_counter + 1;
    fields(field_counter) = {'v_mid'};
    values(:,field_counter) = num2cell(v_mid);
  
% Get n=1 amplitude

  n_equal_1_mode = get_n_equal_1_amplitude(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'n_equal_1_mode'};
    values(:,field_counter) = num2cell(n_equal_1_mode);
  
%}

% Get Te profile width

  Te_HWHM = get_TS_data_cmod(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'Te_width'};
    values(:,field_counter) = num2cell(Te_HWHM);

% Get density and Greenwald fraction

  [n_e, dn_dt, Greenwald_fraction] = get_densities(shot, times_for_db);

    field_counter = field_counter + 1;
    fields(field_counter) = {'n_e'};
    values(:,field_counter) = num2cell(n_e);

    field_counter = field_counter + 1;
    fields(field_counter) = {'dn_dt'};
    values(:,field_counter) = num2cell(dn_dt);

    field_counter = field_counter + 1;
    fields(field_counter) = {'Greenwald_fraction'};
    values(:,field_counter) = num2cell(Greenwald_fraction);

% All the information for this shot has been obtained.  Now prepare to
% enter the information into the disruption_warning table.

  if (field_counter == 0); % If no data is available for any of the fields,
    continue;              % then skip to the next shot
  end;

  fields = fields(1:field_counter);   % Truncate these cell arrays to get
  values = values(:,1:field_counter); % rid of unassigned fields.

% For the time slices that already exist in the database, an SQL 'update'
% operation will be performed.  Otherwise, an SQL 'insert' or 'fastinsert'
% operation will be performed.

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
        if (isnan(cell2mat(values(itime, ifield))));
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
if (strcmp(warning_status.state, 'on'));  % Turn warning messages about calcs
  warning('on', 'MATLAB:interp1:NaNinY'); % w/NaN's back on if it was on
end;                                      % when this routine was called.
