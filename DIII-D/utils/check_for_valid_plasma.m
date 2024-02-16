function [validity, ipmax, duration] = ...
  check_for_valid_plasma(shot, minimum_ip, minimum_duration);

if (~exist('minimum_ip','var')); minimum_ip = 4.e5; end;             % 400 kA
if (~exist('minimum_duration','var')); minimum_duration = 0.1; end;  % 0.1 s

% Check to see if the plasma current signal exists in the D3D tree.  If
% there is no ip data, or no D3D tree for this shot, then the shot is not
% valid.  Exit routine.

mdsconnect('atlas.gat.com'); % connect to mdsserver at DIII-D

[shotopened, status] = mdsopen('d3d', shot);
if (mod(status,2) == 0);  % If d3d tree cannot be opened, the shot is
  validity = 0;           % not valid.  Exit routine.
  ipmax = NaN;
  duration = NaN;
  return;
end;

% Read in Ip data (in amperes)

[ip, status] = mdsvalue(['ptdata("ip", '  int2str(shot)  ')']);
if (mod(status,2) == 0);  % If there is no plasma current data, the shot is
  validity = 0;           % not valid.  Exit routine
  ipmax = NaN;
  duration = NaN;
  return;
end;

ip_time = mdsvalue(['dim_of(ptdata("ip", '  int2str(shot)  '))']);
ip_time = ip_time/1.e3; % convert from ms to s

% Subtract baseline offset

baseline = mean(ip(1:10));
ip = ip - baseline;

% Call my routine "end_of_current" to determine shot duration and maximum
% plasma current.  Note that the last parameter is the Ip threshold for
% determining the end of the shot.  It is not the same thing as the Ip
% threshold for determining if the shot has sufficient plasma current to be
% valid shot for the disruption warning database.

[duration, ipmax] = end_of_current(ip, ip_time, 100e3);

if ((duration < minimum_duration) || (abs(ipmax) < minimum_ip));
  validity = 0;
else;
  validity = 1;
end;
