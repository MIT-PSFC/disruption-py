function [P_ohm, V_loop] = get_P_ohm(shot, timebase);
%
% This function calculates the ohmic power, P_ohm.  The output vector,
% P_ohm, has the same shape as the input vector, "timebase", i.e. if
% "timebase" is a row vector, P_ohm will be a row vector; if "timebase" is
% a column vector, P_ohm will be a column vector.  Note that Matlab calls
% to mdsvalue return column vectors, so this routine works internally with
% column vectors.
%
% We use the following expression to calculate P_ohm:
%   P_ohm = Ip * V_resistive
%      where: V_resistive = V_loop(surf) - V_inductive
%                         = V_loop(surf) - L * dIp/dt
%      and L = L_internal = mu0 * R0 * li/2
%
%   MFLUX node \analysis::top.mflux:v0 is used for V_loop(surf),
%   EFIT  node \efit_aeqdsk:li         is used for li,
%   Magnetics node \Ip                 is used for Ip
%
% and the time derivative of Ip is calculated using Matlab's "gradient"
% function with some boxcar smoothing (+/- 1 ms, acausal).
%
% If the EFIT data or the MFLUX data or the magnetics Ip data is not
% available, then P_ohm is returned as an array of NaN's (Not-a-Number).
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

% Get edge loop voltage

[shotopened, status] = mdsopen('analysis', shot);
if (mod(status,2) == 0);
  P_ohm = NaN(size(timebase));
  V_loop = NaN(size(timebase))
  return;
end;

[V_loop, status] = mdsvalue('\top.mflux:v0');
[V_loop_time, status] = mdsvalue('dim_of(\top.mflux:v0)');
mdsclose;
if (mod(status,2) == 0 || length(V_loop_time) <= 1);
  P_ohm = NaN(size(timebase));
  V_loop = NaN(size(timebase))
  return;
end;

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
    P_ohm = NaN(size(timebase));
    return;
  end;
end;

% Get internal inductance

[li, status] = mdsvalue('\efit_aeqdsk:li');
[efittime, status] = mdsvalue('dim_of(\efit_aeqdsk:li)');
mdsclose;
if (mod(status,2) == 0 || length(efittime) <= 1);
  P_ohm = NaN(size(timebase));
  return;
end;

inductance = 4*pi*1.e-7 * 0.68 * li/2; % For simplicity, we use R0 = 0.68 m,
                                       % but we could use \efit_aeqdsk:rmagx

% Get plasma current from magnetics measurement

[shotopened, status] = mdsopen('magnetics', shot);
if (mod(status,2)==0);
  P_ohm = NaN(size(timebase));
  return;
end;

ip = mdsvalue('\ip');
[magtime, status] = mdsvalue('dim_of(\ip)');
mdsclose;
if (mod(status,2)==0);
  P_ohm = NaN(size(timebase));
  return;
end;

dipdt = gradient(ip, magtime);
dipdt_smoothed = smooth(dipdt, 11); % use 11-point boxcar smoothing

% Now interpolate all the signals onto the time points specified by the
% input vector "timebase".  All the signals used inside this routine are
% column vectors, since mdsvalue returns column vectors.  But the input
% parameter, 'timebase', may be either a column or row vector.  Therefore,
% in order to continue with our column-based calculations, we need to
% create a copy of 'timebase' that is guaranteed to be a column vector.

if (size(timebase,1) > 1);
  timebase_column = timebase;
else;
  timebase_column = transpose(timebase);
end;

V_loop = interp1(V_loop_time, V_loop, timebase_column, 'linear');
inductance = interp1(efittime, inductance, timebase_column, 'linear');
ip = interp1(magtime, ip, timebase_column, 'linear');
dipdt_smoothed = interp1(magtime, dipdt_smoothed, timebase_column, ...
  'linear');

V_inductive = inductance .* dipdt_smoothed;
V_resistive = V_loop - V_inductive;

P_ohm = ip .* V_resistive;

% All of our calculations have used and produced column vectors, since
% mdsvalue returns column vectors.  If the input parameter, 'timebase', is
% a row vector, then make P_ohm a row vector.

if (size(timebase,2) > 1);
  P_ohm = transpose(P_ohm);
  V_loop = transpose(V_loop);
end;

end
