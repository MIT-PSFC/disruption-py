function [validity, ipmax, duration] = ...
  check_for_valid_plasma(shot, minimum_ip, minimum_duration);

if (~exist('minimum_ip','var')); minimum_ip = 1.e5; end;             % 100 kA
if (~exist('minimum_duration','var')); minimum_duration = 0.5; end;  % 0.5 s

% Check to see if the PCRL01 plasma current signal exists in the
% PCS_EAST tree.  If there is no PCRL01 data, or no PCS_EAST tree, then
% the shot is not valid.  Exit routine.

[shotopened, status] = mdsopen('pcs_east', double(shot));
if (mod(status,2) == 0);  % If PCS_EAST tree cannot be opened, the shot is
  validity = 0;           % not valid.  Exit routine.
  ipmax = NaN;
  duration = NaN;
  return;
end;

[ip, status] = mdsvalue('\pcrl01');  % Read in Ip data (in amperes)
if (mod(status,2) == 0);  % If there is no plasma current data, the shot is
  validity = 0;           % not valid.  Exit routine
  ipmax = NaN;
  duration = NaN;
  return;
end;

[ip_time, status] = mdsvalue('dim_of(\pcrl01)');
if (mod(status,2) == 0);  % If there is no Ip timebase data, the shot is
  validity = 0;           % not valid.  Exit routine
  ipmax = NaN;
  duration = NaN;
  return;
end;

% For shots before year 2014, the PCRL01 timebase needs to be shifted by
% 17.0 ms 

if (shot < 44432);
  ip_time = ip_time - 0.0170;
end;

% High-frequency noise spikes on some shots can cause a problem with the
% time derivative and other computations.  Use a median filter to reduce
% the problem.

% In 2016 there seems to be a Matlab license problem that prevents the
% use of the median filter routine at ASIPP, so I have commented it out:

%ip = medfilt1(ip, 5); % Remove noise spikes with median filter

% Subtract baseline offset

baseindices=find(ip_time <= -5.8); % time before any PF supplies turn on
if (numel(baseindices) > 0);
  baseline=sum(ip(baseindices))/length(baseindices);
  ip = ip - baseline;
end;

% Call my routine "end_of_current" to determine shot duration and maximum
% plasma current

[duration, ipmax] = end_of_current(ip, ip_time, minimum_ip);

if ((duration < minimum_duration) || (ipmax < minimum_ip));
  validity = 0;
else;
  validity = 1;
end;
