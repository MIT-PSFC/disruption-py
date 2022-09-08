% Get list of all disruption shots and times in the C-Mod 2016 campaign

db = set_database('logbook');
result = fetch(db, ['select shot, t_disrupt from disruptions where ' ...
  'shot between 1160101000 and 1161231999 order by shot']);
close(db);
disruption_shotlist = int32(cell2mat(result(:,1)));
t_disrupt = cell2mat(result(:,2));
nshots = length(disruption_shotlist);

ndigits = int2str(floor(log10(nshots)) + 1); % (for formatting in the 
                                             %  upcoming fprintf statement)

% All plasma shots in 2016 have the standard EFIT data at 20 ms intervals
% stored in the C-Mod ANALYSIS tree.  We would like to add to this a burst
% of EFIT data at 1 ms intervals for 20 ms before each disruption.  This
% can be accomplished by running EFIT in "burst" mode (mode 10), with one
% "burst" consisting of the standard times, and another burst consisting of
% the higher frequency, pre-disruption times.  These new EFIT data will be
% written to my EFIT18 tree.

for ishot = 1:nshots;
  shot = disruption_shotlist(ishot);
  fprintf(1,['Processing shot %10i  (%' ndigits 'i/%' ndigits ...
    'i  %6.2f%%)\n'], shot, ishot, nshots, ishot/nshots*100);

% Get the list of standard EFIT times in the C-Mod ANALYSIS tree

  mdsopen('analysis', shot);
  efit_standard_times = mdsvalue('\efit_aeqdsk:time');
  mdsclose;

  dt_slow = round(mean(diff(efit_standard_times))*1.e3)/1.e3; % nearest ms
  tstart_slow = round(efit_standard_times(1)*1.e3)/1.e3; % nearest ms
  ntimes_slow = length(efit_standard_times); % See modification to remove
                                             % any overlap with fast-sampled 
                                             % data (follows in a few lines)

% Add in a burst of higher frequency times before a disruption.

  dt_fast = 0.001;
  tstart_fast = t_disrupt(ishot) - 0.020;
  ntimes_fast = int32((t_disrupt(ishot) - tstart_fast)/dt_fast + 1);

% Remove any overlap between the 'slow' times and the 'fast' pre-disruption
% times.  (This was not done for the original C-Mod dataset, which
% consisted of all shots in 2015.)

  ntimes_slow = length(find(efit_standard_times < tstart_fast - 10e-6));

% Create the EFIT18 tree for this shot.  This must be done before running
% EFIT.

  [~, status] = mdsopen('efit18', shot);
  if (mod(status,2) == 0);
    mdsopen('efit18', -1);
    mdstcl(['create pulse ' int2str(shot)]);
    mdsclose;
  else;
    fprintf(1,'  Note: EFIT18 tree already exists for this shot.\n');
  end;

% Create the script file that EFIT will read its inputs from.  Run EFIT in
% burst mode (mode 10), with two bursts: a slow (standard) one, and a fast
% (pre-disruption) one.  (By specifying a negative number for the number of
% bursts, the C-Mod versions of the EFIT executables will sort the combined
% burst timebases into numerical order so that the full timebase will be
% monotonic.)

  fileID = fopen('efit_script_input.txt', 'w');
  fprintf(fileID, '10\n');
  fprintf(fileID, 'efit18\n');
  fprintf(fileID, '\n');  
  fprintf(fileID, '%i, -2\n', shot);
  fprintf(fileID, '%6.1f, %4.1f, %i\n', tstart_slow*1e3, dt_slow*1.e3, ...
    ntimes_slow);
  fprintf(fileID, '%6.1f, %4.1f, %i\n', tstart_fast*1e3, dt_fast*1.e3, ...
    ntimes_fast);
  fclose(fileID);

  status = unix(['/usr/local/cmod/codes/efit/bin/efitdd ' ...
    '< efit_script_input.txt > efit_script_output.txt']);

  if (status == 0);
    unix('rm -f efit_script_input.txt efit_script_output.txt errfil.out');
  end;
end;

clearvars ans db result disruption_shotlist t_disrupt nshots ndigits ...
  ishot shot status efit_standard_times tstart_slow dt_slow ntimes_slow ...
  tstart_fast dt_fast ntimes_fast fileID;
