load('EAST_NEW_SHOTS_2018_2020.mat');
shotlist = SHOTLIST;
clearvars -except shotlist

mdsconnect('mds.ipp.ac.cn'); % Connect to EAST MDSplus server

nshots = length(shotlist);
ndigits = int2str(floor(log10(nshots)) + 1); % (for formatting)

valid_id = zeros(nshots,1);

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

  if (~valid); 
    continue; 
  else
    valid_id(ishot) = 1;
  end;  % If not a valid shot, skip it


end

valid_plasmas = shotlist(valid_id==1);
save('valid_plasmas_2018-2020.mat', 'valid_plasmas');
