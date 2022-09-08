% Steve Scott's two shots (1101201019 and 1101201020) have the standard
% EFIT data at 20 ms intervals stored in the C-Mod ANALYSIS tree.  The
% first one disrupted, the 2nd one didn't.  I want to rerun EFIT on both of
% these shots at 1 ms intervals for the duration of each shot.  However,
% there is a limit of 301 reconstructions per shot.  So I will use 1 ms
% intervals for the 200 ms period before the disruption, and I will use the
% same 200 ms high speed sampling window for the non-disruptive shot as
% well.  The rest of the shot will use the original 20 ms EFIT timing.  The
% new EFIT data will be written to my EFIT18 trees.

shotlist = [1101201019, 1101201020];
t_disrupt = [1.1195, 1.3095]; % The 2nd shot didn't disrupt, but we want to
nshots = length(shotlist);    % analyze it as if it did.

for ishot = 1:nshots;
  shot = shotlist(ishot);
  fprintf(1,['Processing shot %10i  (%i/%i ' ...
    '  %6.2f%%)\n'], shot, ishot, nshots, ishot/nshots*100);

% Get the list of standard EFIT times in the C-Mod ANALYSIS tree

  mdsopen('analysis', shot);
  efit_standard_times = mdsvalue('\efit_aeqdsk:time');
  mdsclose;

% First, do the times before the window of high sampling rate

  indices_standard_times_before = ...
    find(efit_standard_times < t_disrupt(ishot) - 0.200 - 10e-6);
  efit_standard_times_before = ...
    efit_standard_times(indices_standard_times_before);

  dt_before = ...
    round(mean(diff(efit_standard_times_before))*1.e3)/1.e3; % nearest ms
  tstart_before = ...
    round(efit_standard_times_before(1)*1.e3)/1.e3; % nearest ms
  ntimes_before = length(efit_standard_times_before);

% Add in a burst of higher frequency times before a disruption.

  dt_fast = 0.001;
  tstart_fast = t_disrupt(ishot) - 0.200;
  ntimes_fast = int32((t_disrupt(ishot) - tstart_fast)/dt_fast + 1);

% Finally, do the times after the window of high sampling rate.  For the
% shot that disrupts, there shouldn't be any times after the fast window.

  indices_standard_times_after = ...
    find(efit_standard_times > t_disrupt(ishot) + 10e-6);
  efit_standard_times_after = ...
    efit_standard_times(indices_standard_times_after);

  dt_after = ...
    round(mean(diff(efit_standard_times_after))*1.e3)/1.e3; % nearest ms
  tstart_after = ...
    round(efit_standard_times_after(1)*1.e3)/1.e3; % nearest ms
  ntimes_after = length(efit_standard_times_after);

  ntimes = ntimes_before + ntimes_fast + ntimes_after;
  if (ntimes > 301);
    ntimes_after = ntimes_after - (ntimes - 301);
  end;

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
  if (ntimes_after > 1);
    fprintf(fileID, '%i, -3\n', shot);
  else;
    fprintf(fileID, '%i, -2\n', shot);
  end;
  fprintf(fileID, '%6.1f, %4.1f, %i\n', tstart_before*1e3, ...
    dt_before*1.e3, ntimes_before);
  fprintf(fileID, '%6.1f, %4.1f, %i\n', tstart_fast*1e3, ...
    dt_fast*1.e3, ntimes_fast);
  if (ntimes_after > 1);
   fprintf(fileID, '%6.1f, %4.1f, %i\n', tstart_after*1e3, ...
      dt_after*1.e3, ntimes_after);
  end;
  fclose(fileID);

  status = unix(['/usr/local/cmod/codes/efit/bin/efitdd ' ...
    '< efit_script_input.txt > efit_script_output.txt']);

  if (status == 0);
    unix('rm -f efit_script_input.txt efit_script_output.txt errfil.out');
  end;

end;
