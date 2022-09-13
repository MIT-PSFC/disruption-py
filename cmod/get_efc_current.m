function Iefc = get_efc_current(shot, timebase);

% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Output:
%   Iefc = current from EFC chopper power supply vs time [A]
%
% Author: Robert Granetz  June 2020
% Update history:

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
else
  timebase_column = transpose(timebase);
end;

% Initialize all output arrays to NaN (Not-a-Number).

Iefc = NaN(size(timebase));

% Read in the EFC current from the engineering tree. 

[~,status] = mdsopen('engineering',shot);
if (mod(status,2) == 1)
  Iefc = mdsvalue('\efc:u_bus_r_cur');
  [time, status] = mdsvalue('dim_of(\efc:u_bus_r_cur)');
  if (mod(status,2) ~= 1 || length(Iefc) ~= length(time));
    return;
  end

  if (mod(status,2) == 1)
      Iefc = interp1(time, Iefc, timebase_column, 'linear');
  end  
end
mdsclose; 

% The output signals are currently all column vectors.  However, we desire
% to have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase,2) > 1);
  Iefc = transpose(Iefc);
end;

end
