function [Z_cur, vZ, Z_times_vZ] = get_Z(shot, timebase)

% This script reads in Z_cur from EFIT, then calculates the vertical
% velocity, vZ, and the product, Z_cur * vZ.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   Z_cur = EFIT's reconstruction of the vertical position [m]
%   vZ = vertical velocity [m/s]
%   Z_times_vZ = Z_cur1*vZ [m^2/s]
%
% Author: Alex Tinguely 2015-09-09
% Modified by Robert Granetz on 2016/01/28

% The input array, "timebase", can be either a row vector (1 x n) or a
% column vector (n x 1).  We desire to have the output arrays match the
% shape of "timebase".  However, Matlab calls to the routine "mdsvalue"
% return column vectors for 1-D signals, so it is simpler to work with
% column vectors within this routine, and then, if necessary, convert the
% outputs to row vectors just before exiting this routine.  So the first
% step is to create a copy of "timebase" that is guaranteed to be a column
% vector.

if (size(timebase,1) > 1);
  timebase_column = timebase;
else;
  timebase_column = transpose(timebase);
end;

% Initialize all output arrays to NaN (Not-a-Number), with the same shape
% as "timebase".

Z_cur = NaN(size(timebase));
vZ = NaN(size(timebase));
Z_times_vZ = NaN(size(timebase));

% EFIT data on the standard timebase are in the ANALYSIS tree.  For shots
% that disrupted, these data are superceded by EFIT data in the EFIT18
% tree, which includes data at a high sampling rate prior to the disruption
% time, in addition to the standard times.  This routine does not know, a
% priori, which tree it should open.  Therefore, I will first try to open
% the EFIT18 tree, and if this fails, then I'll open the ANALYSIS tree.

[shotopened, status] = mdsopen('efit18', shot);
if (mod(status,2) == 0);
  [shotopened, status] = mdsopen('analysis', shot);
  if (mod(status,2)==0);
    return;
  end;
end;

% Read in EFIT's reconstruction of the vertical position, and its timebase

[Z_cur, status] = mdsvalue('\efit_aeqdsk:zcur');
[efittime, status] = mdsvalue('dim_of(\efit_aeqdsk:zcur)');
mdsclose;

if (mod(status, 2) == 0); % if an error occurs
  return;
end;

vZ = gradient(Z_cur, efittime); % vZ = d(Z_cur)/dt
Z_times_vZ = Z_cur .* vZ; 

% Interpolate signals onto the specified timebase

Z_cur = interp1(efittime, Z_cur, timebase_column, 'linear');
vZ = interp1(efittime, vZ, timebase_column, 'linear');
Z_times_vZ = interp1(efittime, Z_times_vZ, timebase_column, 'linear');

% If "timebase" is a row vector, then convert all the outputs back to row
% vectors

if (size(timebase,2) > 1);
  Z_cur = transpose(Z_cur);
  vZ = transpose(vZ);
  Z_times_vZ = transpose(Z_times_vZ);
end;

end