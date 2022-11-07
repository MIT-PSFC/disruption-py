function [v_0, v_mid] = get_rotation_velocity(shot, times_for_db)

% This routine obtains toroidal rotation velocities (on-axis and 
% mid-radius). 
% 
% Inputs:
%   shot = shot number
%   time_for_db = array of desired time values?
%
% Outputs:
%   v_0 = toroidal rotation velocity on-axis
%   v_mid = toroidal rotation velocity at mid-radius
%
% Authors: Kevin Montes and Robert Granetz   Jan 2017

% The input array, "times_for_db", can be either a row vector (1 x n) or a
% column vector (n x 1).  We desire to have the output arrays match the
% shape of "times_for_db".  However, Matlab calls to the routine "mdsvalue"
% return column vectors for 1-D signals, so it is simpler to work with
% column vectors within this routine, and then, if necessary, convert the
% outputs to row vectors just before exiting this routine.  So the first
% step is to create a copy of "times_for_db" that is guaranteed to be a
% column vector.

if (size(times_for_db,1) > 1);
  times = times_for_db;
else
  times = transpose(times_for_db);
end;

% Initialize all output arrays to NaN (Not-a-Number) column vectors

v_0   = NaN(size(times));
v_mid = NaN(size(times));

% Check to see if shot was done on a day where there was a locked
% mode HIREX calibration by cross checking with list of calibrated 
% runs. If not calibrated, return NaN outputs.

list = importdata( ...
  '/home/montes/Disruption_warning_code/lock_mode_calib_shots.txt');
run = floor(shot/1000);
if ~any(list == run)
    return
end

% Next, check that (1) the spectroscopy tree opens, (2) intensity
% and velocity data are stored, and (3) the argon intensity pulse 
% meets a minimum count and duration threshold. If (1)-(3) are
% satisfied, linearly interpolate from HIREX timebase to input times
% and assign velocity values to v_0. If any of (1)-(3) are not 
% satisfied, return all v_0 values as NaN. 

node = '.hirex_sr.analysis.a';

T = 0.2; % Min. intensity pulse duration (secs) b/w 10^3,10^4 counts
[~, status] = mdsopen('spectroscopy',shot);
if (mod(status,2) == 1); % (1)
    intensity = mdsvalue([node ':int']);
    [~, status] = mdsvalue(['dim_of(' node ':int)']);
    if (mod(status,2) == 1); % (2)
        vel = mdsvalue([node ':vel']);
        [hirextime, status] = mdsvalue(['dim_of(' node ':vel)']);
        if (mod(status,2) == 1); % (2)
            indices = find(intensity>1000 & intensity <10000);
            goodtime = length(indices)*.02;
            if goodtime >= T % (3)
                v_0 = interp1(hirextime, vel, times, 'linear');
            end               
        else
            return
        end
    else
        return
    end
    mdsclose;
else
    return
end

% Exclude any velocities greater than 200 km/s and convert all other
% velocities to m/s. 

for i=1:size(v_0,1)
    if abs(v_0(i)) > 200 % If I find stricter condition to exclude,put here
        v_0(i) = NaN;
    end
end

v_0 = v_0 * 1000;
