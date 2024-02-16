function [p_RAD, p_NBI, p_OHM, p_ECH, rad_fraction, p_input, v_loop] = ...
    get_power_d3d(shot, timebase);

% This function gets the input heating powers -- ohmic (p_OHM), electron
% cyclotron (p_ECH), and neutral beam (p_NBI) -- as well as the radiated
% output power (p_RAD).  If any of the auxiliary heating powers are not
% available (e.g. there was no ECH), then this function returns an array of
% zeros for them.  The ohmic heating power is obtained by calling a
% separate function, get_P_ohm_d3d.m.  If either the ohmic heating power or
% the radiated power is unavailable, then arrays of NaN (Not-a-Number) are
% returned for them, and for the radiated power fraction.
%
% Inputs:
%   shot = shot number
%   timebase = times to calculate the powers 
% Outputs:
%   p_RAD  = radiated power [W]
%   p_NBI  = neutral beam power [W]
%   p_OHM  = ohmic power [W]
%   p_ECH  = electron cyclotron power [W]
%   p_input = total input power = p_OHM + p_NBI + p_ECH [W]
%   rad_fraction = p_RAD/(p_OHM + p_NBI + p_ECH), ratio of radiated output
%       power to input heating power 
%
% NOTE: In DIII-D, there sometimes can be ICH (ion cyclotron heating)
% power, but this is rare and usually negligible compared to the neutral
% beam injection.
%
% Author: Alex Tinguely 15-09-16
% Some modifications by Robert Granetz 2015/10/11
% Updated by Alex Tinguely for DIII-D 2016-02-29 (Leap Day!)
% Updated by Cristina Rea for DIII-D 2016-09-19 (implemented automated
%   selection of efit trees)
% Some modifications by Robert Granetz 2016/11/11
% Updated by Kevin Montes 2020-05-05 (fixed error retrieving ECH power)

% We want the output vectors to have the same shape as the input vector,
% "timebase", i.e. if "timebase" is a row vector, the output vectors will
% be row vectors; if "timebase" is a column vector, the output vectors will
% be column vectors.  However, since the Matlab version of 'mdsvalue'
% returns column vectors for 1-D signals, this routine works internally
% with column vectors.  Just before this routine returns, it transforms the
% column vectors to match the shape of the "timebase" input vector, if
% necessary.  So our first step in this routine is to create a column
% vector copy of the input parameter, "timebase".

if (size(timebase,1) > 1);
  timebase_column = timebase;
else;
  timebase_column = transpose(timebase);
end;

% ------------------------------------
% Get neutral beam injected power

[shotopened, status] = mdsopen('d3d', shot);

if (mod(status, 2)==1);
  p_nbi = mdsvalue('\d3d::top.nb:pinj'); % [kW]
  p_nbi = p_nbi * 1.e3; % convert to watts
  [t_nbi, status] = mdsvalue('dim_of(\d3d::top.nb:pinj)'); % [ms]
  t_nbi = t_nbi/1.e3; % convert to seconds
  if (mod(status,2)==1 && length(t_nbi) > 2);
    p_NBI = interp1(t_nbi, p_nbi, timebase_column, 'linear', 0.);
  else;
    p_NBI = zeros(size(timebase_column));
  end;
  mdsclose;
else;
  p_NBI = zeros(size(timebase_column));
end;

% ------------------------------------
% Get ECH power.  It's point data, so it's not stored in an MDSplus tree

[shotopened, status] = mdsopen('rf', shot);

if (mod(status, 2)==1);
  p_ech = mdsvalue('\top.ech.total:echpwrc'); % [W]
  [t_ech, status] = mdsvalue('dim_of(\top.ech.total:echpwrc)'); % [ms]
  t_ech = t_ech/1.e3; % convert to seconds
  if (mod(status,2)==1 && length(t_ech) > 2);
	% Sometimes, t_ech has an extra "0" value tacked on to the end. This must
	% be removed before the interpolation.
	if t_ech(end)==0
		t_ech = t_ech(1:end-1); p_ech = p_ech(1:end-1);
  	end
    p_ECH = interp1(t_ech, p_ech, timebase_column, 'linear', 0.);
  else;
    p_ECH = zeros(size(timebase_column));
  end;
  mdsclose;
else;
  p_ECH = zeros(size(timebase_column));
end;

%{
%% Old ECH code using echpwr pointname (before 05/2020 update)
%-------------------------------------------------------------
  p_ech = mdsvalue(['ptdata("echpwr", ' num2str(shot) ')']); % [W]
  [t_ech, status] = ...
    mdsvalue(['dim_of(ptdata("echpwr", ' num2str(shot) '))']); % [ms]
  t_ech = t_ech/1.e3; % convert to seconds
    
  if (mod(status,2)==1 && length(t_ech) > 1);

% Sometimes, the length of t_ech is not the same as p_ech, and this seems
% to be due to an extra "0" value tacked onto the end of t_ech (right after
% the 10.0 s value), so in this case we remove the last value.

    if length(t_ech) ~= length(p_ech);
       t_ech = t_ech(1:end-1); 
    end;
    p_ECH = interp1(t_ech, p_ech, timebase_column, 'linear', 0.);
  else;
    p_ECH = zeros(size(timebase_column));
  end;
%}

% ------------------------------------
% Get ohmic power and loop voltage

[p_OHM, v_loop] = get_P_ohm_d3d(shot, timebase); %[W, V]

% ------------------------------------
% Radiated power
%
% We had planned to use the standard signal '\bolom::prad_tot' for this
% parameter.  However, the processing involved in calculating \prad_tot
% from the arrays of bolometry channels involves non-causal filtering with
% a 50 ms window.  This is not acceptable for our purposes.  Tony Leonard
% provided us with the two IDL routines that are used to do the automatic
% processing that generates the \prad_tot signal in the tree (getbolo.pro
% and powers.pro).  I converted them into Matlab routines, and modified the
% analysis so that the smoothing is causal, and uses a shorter window.

smoothing_window = 0.010; % use 10 ms causal smoothing window
a_structure = getbolo_new(shot, smoothing_window*1.e3);

% Sometimes the bolo data is garbage.  Check the 'ier' subfields to
% determine this

ier = 0;
for ichan = 1:48;
  if a_structure.chan(ichan).ier == 1;
    ier = 1;
  end;
end;
if ier == 1;
  p_RAD = NaN(size(timebase_column));
else;
  b_structure = powers_new(a_structure);
  p_RAD = b_structure.pwrmix; % watts
  p_RAD_timebase = a_structure.rawtime; % seconds

% Make sure that p_RAD and its timebase are column vectors
  if size(p_RAD, 2) > 1;
    p_RAD = transpose(p_RAD);
    p_RAD_timebase = transpose(p_RAD_timebase);
  end;

% Interpolate onto the desired timebase
  p_RAD = interp1(p_RAD_timebase, p_RAD, timebase_column, 'linear');
end;

% ------------------------------------
% Get rid of negative values of p_RAD, p_NBI, and p_ECH (unphysical).

p_RAD(p_RAD < 0) = 0;
p_NBI(p_NBI < 0) = 0;
p_ECH(p_ECH < 0) = 0;

% If "timebase" is a row vector, then transform all the outputs to row
% vectors.  Otherwise, leave them as column vectors.  Note that p_OHM and
% v_loop are already correct, since that is done in the get_P_ohm_d3d
% routine.

if (size(timebase,2) > 1);
  p_NBI = transpose(p_NBI);
  p_ECH = transpose(p_ECH);
  p_RAD = transpose(p_RAD);
end;

p_input = p_OHM + p_NBI + p_ECH; % Total input power [W]
rad_fraction = p_RAD ./ p_input; % Ratio of radiated power to input power

rad_fraction(isinf(rad_fraction)) = NaN; % Change any infinite values
                                         % to Not-a-Number
end
