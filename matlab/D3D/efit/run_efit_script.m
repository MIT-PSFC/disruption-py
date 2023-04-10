cd('/fusion/projects/disruption_warning/software/matlab_programs/efit_runs/');
load('d3d_shot_duration_disrupt.mat'); % shotlist, end_of_shot, disrupt_time

delta_t_no_disrupt = 0.025;  % Do EFIT's every 25 ms on all plasmas.
delta_t_disrupt    = 0.002;  % For disruptions, also do EFIT's every 2 ms,
disrupt_window     = 0.100;  % starting 0.100 seconds before the disruption.

%for ishot = 1:length(shotlist);
%for ishot = 2054:length(shotlist);
for ishot = 2853:2999;
  shot = shotlist(ishot);
  if (isnan(disrupt_time(ishot)));
    fprintf(1, 'Processing shot %i  [\43%4i]\n', shot, ishot);
  else;
    fprintf(1, 'Processing shot %i  [\43%4i] (disrupted)\n', shot, ishot);
  end;

% Create subdirectory for this shot

  shotdir = ['shot' num2str(shot,'%i')];
  mkdir(shotdir);

% Create EFIT snapfile for this shot

  create_efit_snapfile(shot);

% Create EFIT command input file for this shot.  This file contains the
% specification of the time slices to do.

  cd(shotdir);
  ntimes = round((end_of_shot(ishot)- 0.100)/delta_t_no_disrupt) + 1;
  fileid = fopen('efit_infile.txt', 'w');
  fprintf(fileid, '7\n');
  fprintf(fileid, 'disrupt\n');
  fprintf(fileid, '%i, %i, %i, %i\n', shot, round(0.100*1.e3), ...
    round(delta_t_no_disrupt*1.e3), min([ntimes, 400]));
  fprintf(fileid, '\n');
  fclose(fileid);

% Run EFIT, reading its input from 'efit_infile.txt'

  [status_efit, cmdout_efit] = unix('efitd90 65 65 < efit_infile.txt');
%{
  if (~status_efit);
    fprintf(1,' EFIT returned error to shell on shot %i\n', shot);
    fprintf(1, '%s\n', cmdout_efit);
  end;
%}

% If this shot disrupted, run EFIT with higher sampling rate near the
% disruption time.  First I have to create a new 'efit_infile.txt' file
% that specifies the disruption-relevant time slice parameters.

  if (~isnan(disrupt_time(ishot)));
    tdis = disrupt_time(ishot);
    ntimes = round(disrupt_window/delta_t_disrupt) + 1;
    ntimes = min([ntimes, 400]);
    fileid = fopen('efit_infile.txt', 'w');
    fprintf(fileid, '7\n');
    fprintf(fileid, 'disrupt\n');
    fprintf(fileid, '%i, %i, %i, %i\n', shot, ...
      round((tdis-disrupt_window)*1.e3), round(delta_t_disrupt*1.e3), ...
      min([ntimes, 400]));
    fprintf(fileid, '\n');
    fclose(fileid);

    [status_efit, cmdout_efit] = unix('efitd90 65 65 < efit_infile.txt');
%{
    if (~status_efit);
      fprintf(1,' EFIT returned error to shell on shot %i\n', shot);
      fprintf(1, '%s\n', cmdout_efit);
    end;
%}
  end;

% The EFIT runs should have created eqdsk files for each requested time slice.

% Now run EFITLOADER to collect all the eqdsk time slices, put then in
% temporal order, and write them into a new EFIT tree.

% First, I have to create the "efitloader.info" file for the current shot.

  fileid = fopen('efitloader.info', 'w');
  fprintf(fileid, '4\n');
  fprintf(fileid, ['FILESPEC="/fusion/projects/disruption_warning/' ...
                   'software/matlab_programs/efit_runs/shot%i"\n'], shot);
  fprintf(fileid, 'NAMELIST="efit_snap.dat_disrupt"\n');
  fprintf(fileid, 'COMMENTS="for disruption warning database"\n');
  fprintf(fileid, 'RUNTAG="DIS"\n');
  fclose(fileid);

% Okay, now run EFITLOADER

  cmd = ['efitloader -mode PASTED -noarchive -nodialog -full ' ...
         '-f "/fusion/projects/disruption_warning/software/matlab_programs/' ...
         'efit_runs/shot' num2str(shot) '/efitloader.info"'];
  [status_efitloader, cmdout_efitloader] = unix(cmd);

  cd('/fusion/projects/disruption_warning/software/matlab_programs/efit_runs/');

% Now that the EFIT tree has been created, I can delete all the EFIT files
% and the shot directory.

  [status] = unix(['rm -f ' shotdir '/*']);
  rmdir(shotdir, 's');

end;
