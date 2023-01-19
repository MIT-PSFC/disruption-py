function V_loop = get_vloop_d3d(shot, timebase);

% This function calculates the ohmic power, P_ohm, and also outputs the
% loop voltage, V_loop.
%
% We use the following expression to calculate P_ohm:
%   P_ohm = Ip * V_resistive
%      where: V_resistive = V_loop(surf) - V_inductive
%                         = V_loop(surf) - L * dIp/dt
%      and L = L_internal = mu0 * R0 * li/2
%
% For DIII-D:
%   ptdata("vloopb") in d3d is used for V_loop (measured outside the vessel),
%   EFIT node \efit_a_eqdsk:li is used for li,
%   ptdata("ip") in d3d is used for Ip
%
% and the time derivative of Ip is calculated using Matlab's "gradient"
% function, with some boxcar smoothing (10 ms, CAUSAL).
%
% If the EFIT data or the V_loop data or the Ip data is not available, then
% P_ohm is returned as an array of NaN's (Not-a-Number).
%
% Inputs:
%   timebase = times to calculate P_ohm and V_loop at 
% Outputs:
%   P_ohm = ohmic power [W]
%   V_loop = loop voltage [V]
%
% Author: Alex Tinguely 15-09-16
% Some modifications by Robert Granetz 2015/10/10
%                                      2016/01/28
%                                      2016/02/09 (changed Vloop to MFLUX V0)
% Updates for DIII-D by Alex Tinguely on 2016-02-26
% Updates for DIII-D by Cristina Rea on 2016-09-19 (implemented automated
%   selection of efit trees)
% Some modifications by Robert Granetz 2016/11/11 (use ptdata("vloopb") for
%   the loop voltage; use 10 ms causal filter for smoothing dip/dt

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

% Get edge loop voltage and smooth it a bit with a median filter

[V_loop, status] = mdsvalue(['ptdata("vloopb", ' num2str(shot) ')']); % Volts
[V_loop_time, status] = ...
   mdsvalue(['dim_of(ptdata("vloopb", ' num2str(shot) '))']); % ms
V_loop_time = V_loop_time/1000; % convert to seconds
% median filter; signal distortion increases by increasing the order number
% need to use odd-number filtering to avoid “forward-shifting”;
% vloopb is sampled at 2kHz, we choose a window of 5ms for the median
% filter
order = 11;
V_loop = medfilt1(V_loop, order); 

if (mod(status,2) == 0 || length(V_loop_time) <= 1);
  P_ohm = NaN(size(timebase));
  V_loop = NaN(size(timebase));
  return;
else;
  V_loop = interp1(V_loop_time, V_loop, timebase_column, 'linear');
end;

% All of our calculations have used and produced column vectors, since
% mdsvalue returns column vectors.  If the input parameter, 'timebase', is
% a row vector, then transpose the outputs into row vectors in order to
% match the shape of 'timebase'.

if (size(timebase,2) > 1);
  V_loop = transpose(V_loop);
end;

end
