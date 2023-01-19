function create_d3d_shot_duration_disrupt_datafile;

%  2022-04-25 CR: this routine returns the .mat file needed
%                 to run the main disruption_warning_database.m
%                 AND to run modified EFITs
%                 Modified routine to read in 'shot' from SQL db
%  INPUT: shot is the last shot # currently in the db


% add path to find end_of_current.m
addpath('/fusion/projects/disruption_warning/software/matlab_programs/');

db = set_database('d3drdb');

% added CR
last_shot = cell2mat(fetch(db,['select shot from disruption_warning order by shot']));
shot = max(unique(last_shot));
%

shots= cell2mat(fetch(db,['select shot from summaries where ' ...
  'shot >= ' int2str(shot) ' and pulse_length > 0.1 and abs(ip) > 400e3 order by shot']));

results = fetch(db,['select shot, t_disrupt from disruptions where ' ...
  'shot >= ' int2str(shot) ' order by shot']);
close(db);

disrupted_shots = cell2mat(results(:,1));
t_disrupt = cell2mat(results(:,2));
clearvars results db;

mdsconnect('atlas.gat.com');

shotlist = unique(shots);
nshots = length(shotlist);
ndigits = int2str(floor(log10(nshots)) + 1); % (for formatting in the 
                                             %  upcoming fprintf statement)
end_of_shot = NaN(size(shotlist));
disrupt_time = NaN(size(shotlist));

for ishot = 1:length(shotlist);

  this_shot = shotlist(ishot);
  fprintf(1,['Processing shot %7i  (%' ndigits 'i/%' ndigits ...
    'i  %6.2f%%)\n'], this_shot, ishot, nshots, ishot/nshots*100);

  [shotopened, status]=mdsopen('d3d', this_shot);
  [iptime, status] = mdsvalue(['dim_of(ptdata("ip", ' num2str(this_shot) '))']);
  if (mod(status,2) == 1);
    iptime = iptime/1.e3; % convert ms to s
    ip = mdsvalue(['ptdata("ip", ' num2str(this_shot) ')']);
    mdsclose;
  else;
    fprintf(1,'  Unable to open D3D tree for shot%7i\n', this_shot);
  end;
  [times, Ipmax] = end_of_current(ip,iptime,1.e5);
  end_of_shot(ishot) = max(times);

  [disrupted,indx] = ismember(this_shot, disrupted_shots);
  if (disrupted);
    disrupt_time(ishot) = t_disrupt(indx);
  end; 
end;

save('d3d_shot_duration_disrupt.mat', 'shotlist', 'end_of_shot', ...
  'disrupt_time');
mdsdisconnect;
end
