function [P_ohm, V_loop] = get_P_ohm_d3d(shot, timebase);

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

% Get plasma current

ip = mdsvalue(['ptdata("ip", ' num2str(shot) ')']);
[iptime, status] = mdsvalue(['dim_of(ptdata("ip", ' num2str(shot) '))']);

if (mod(status,2) == 1);
  iptime = iptime/1.e3; % convert ms to s

%{
  dipdt = gradient(ip, iptime);

%  Filter dip/dt using CAUSAL moving average filter (i.e. boxcar).  The Ip
%  signal on DIII-D is a lot noisier than on C-Mod, so I need to use a
%  wider smoothing window, but at least it's causal, and the ohmic power is
%  usually a small fraction of the total heating power anyway during the
%  flattop.

  boxcar_width = 0.010; % use 10 ms wide boxcar window for causal Ip smoothing
  dt = iptime(2)-iptime(1);
  windowsize = round(boxcar_width / dt);
  smoothing_kernal = (1/windowsize) * ones(1, windowsize);
  dipdt_smoothed = filter(smoothing_kernal, 1, dipdt);
%}

% Alessandro Pau (JET & AUG) has given Cristina a robust routine that
% performs time differentiation with smoothing, while preserving causality.
% It can be useful for differentiating numerous signals such as Ip, Vloop,
% etc.  It is called 'GSASTD', and Cristina has put it in our disruption
% warning project directory.  We will use this routine in place of Matlab's
% 'gradient' and smoothing/filtering routines.
%      Processed = GSASTD(x,y,DerivativeMode,width,type,ends,SlewRate)
% we choose a 20-points width for GSASTD: this means a 10ms window for Ip
% smoothing
  width = 20;
  dipdt_smoothed = GSASTD(iptime, ip, 1, width, 3, 1, 0);
else;
  P_ohm = NaN(size(timebase));
  if (size(timebase,2) > 1); V_loop = transpose(V_loop); end;
  return;
end;

% We have rerun EFIT for most DIII-D discharges in the 2015 campaign at our
% desired sampling times, in order to avoid interpolating EFIT signals, and
% also to use minimal non-causal filtering of the magnetics data.  These
% EFIT tree numbers get assigned by an archiving system, and vary from
% shot-to-shot.  We want to select the highest EFIT tree number that was
% created by 'granetzr' and has the run label 'DIS'.  If no such EFIT tree
% exists, then use the default post-shot tree, EFIT01.

% First, turn off annoying java warning messages that tend to be generated
% when using the SQL java driver

warning_status = warning;
warning('off')

efittrees = select_efit_trees(shot,'','DIS');

% Return java warning message status to original setting
warning(warning_status);

if ~isempty(efittrees);
  tree=char(efittrees(end,:));
else;
  fprintf(1, 'No disruption EFIT tree for this shot.  Using EFIT01.\n');
  tree = 'EFIT01';
end;

% Open EFIT tree and get internal inductance

[shotopened, status] = mdsopen(tree, shot);
if (mod(status,2) == 1);
  [li, status] = mdsvalue('\efit_a_eqdsk:li');
  [efittime, status] = mdsvalue('dim_of(\efit_a_eqdsk:li)'); % in ms
  efittime = efittime/1000; % convert to seconds
  chisq = mdsvalue('\efit_a_eqdsk:chisq'); % Use chisq to determine which time
                                         % slices are invalid.
  mdsclose;

  % EFIT reconstructions are sometimes invalid, particularly when very close
% to a disruption.  There are a number of EFIT parameters that can indicate
% invalid reconstructions, such as 'terror' and 'chisq'.  Here we use
% 'chisq' to determine which time slices should be excluded from our
% disruption warning database.

  invalid_indices = find(chisq > 50);
  li(invalid_indices) = NaN;

  if (mod(status,2) == 0 || length(efittime) <= 1);
    P_ohm = NaN(size(timebase));
    if (size(timebase,2) > 1); V_loop = transpose(V_loop); end;
    return;
  end;
else;
  P_ohm = NaN(size(timebase));
  if (size(timebase,2) > 1); V_loop = transpose(V_loop); end;
  return;
end;

R0 = 1.67; % For simplicity, use fixed R0 = 1.67 m for DIII-D major radius
inductance = 4*pi*1.e-7 * R0 * li/2; % henries

inductance = interp1(efittime, inductance, timebase_column, 'linear');
ip = interp1(iptime, ip, timebase_column, 'linear');
dipdt_smoothed = interp1(iptime, dipdt_smoothed, timebase_column,'linear');

V_inductive = inductance .* dipdt_smoothed;
V_resistive = V_loop - V_inductive;

P_ohm = ip .* V_resistive;

% All of our calculations have used and produced column vectors, since
% mdsvalue returns column vectors.  If the input parameter, 'timebase', is
% a row vector, then transpose the outputs into row vectors in order to
% match the shape of 'timebase'.

if (size(timebase,2) > 1);
  P_ohm = transpose(P_ohm);
  V_loop = transpose(V_loop);
end;

end
