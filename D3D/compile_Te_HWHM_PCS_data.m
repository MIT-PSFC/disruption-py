
% Revision/history:
%  2017/08/22 - RSG; Modify this program to only compile Thomson profile
%                    widths using the actual PCS data, rather than the
%                    blessed data in the MDS trees.  Compile these data
%                    for a set of 392 shots selected by Cristina Rea, and
%                    store them in a temporary file.

% Define path on the Iris cluster to some required Matlab routines

oldpath = addpath('/fusion/projects/disruption_warning/matlab_programs');

% Turn off annoying warning messages about calculations with NaN values in
% the Matlab INTERP1 interpolation routine.

warning_status = warning;
warning('off')

% Connect to the SQL database that contains the 'disruption_warning' table

db = set_database('d3drdb'); % Connect to the SQL database that contains
                             % the 'disruptions' table

% Get list of shots selected by Cristina Rea

shotlist=dlmread(['/fusion/projects/disruption_warning/matlab_programs/' ...
  'shotlist_rea_blessed.txt']);
shotlist = int32(shotlist);

% Create array to hold the TS profile width data

Te_HWHM_PCS = NaN(0);

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

% Get list of time slices in the database for this shot

  result = fetch(db, ['select time from disruption_warning ' ...
    'where shot = ' int2str(shot) ' order by time']);
  timebase = cell2mat(result(:,1));
  ntimes = length(timebase);

% Get Te profile width (from Thomson scattering data)

  [Te_HWHM_PCS_this_shot] = get_TS_data_PCS( shot, timebase);

% Add this shot's data to the array

  Te_HWHM_PCS = [Te_HWHM_PCS; Te_HWHM_PCS_this_shot];

end;

save('Te_HWHM_PCS.mat', 'Te_HWHM_PCS');

% Okay, we have processed all the shots in the list.  We just need to clean
% up a few things, and then we're done.

path(oldpath);
close(db);

%{
if (strcmp(warning_status.state, 'on'));  % Turn warning messages about calcs
  warning('on', 'MATLAB:interp1:NaNinY'); % w/NaN's back on if it was on
end;                                      % when this routine was called.
%}

warning(warning_status);
