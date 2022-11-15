function [beta_N, beta_p, dbeta_p_dt, kappa, upper_gap, lower_gap, ...
  li, dli_dt, q0, qstar, q95, V_loop_efit, Wdiam, dWdiam_dt] = ...
  get_EFIT_parameters(shot, timebase);

% This function obtains the EFIT parameters, and also calculates the time
% derivatives of a few of them.
%
% Inputs:
%   shot = shot number
%   timebase = array of times at which to evaluate the parameters
%
% Outputs:
%   beta_N = normalized beta
%   beta_p = beta_poloidal
%   dbeta_p_dt = time derivative of beta_p [/s]
%   kappa = elongation
%   upper_gap = upper gap [m]
%   lower_gap = lower gap [m]
%   li = internal inductance
%   dli_dt = time derivative of li [/s]
%   q0 = safety factor at center
%   qstar = q* = cylindrical safety factor
%   q95 = edge safety factor
%   V_loop_efit = loop voltage [V]
%   Wdiam = diamagnetic energy, [J]
%   dWdiam_dt = time derivative of Wdiam, [J/s]
%
% Author: Alex Tinguely 2015-10-07 (Dan Brunner's 30th Birthday!)
% Some modifications by Robert Granetz 2015/10/20
% Updates: DIII-D by Alex Tinguely 2016-02-26
% Updates: DIII-D EFIT parameters changes by Cristina Rea 2016-08-11

% Initialize all outputs to NaN (Not-a-Number).  These will be returned by
% this routine if it is unable to open the MDSplus tree or read out the
% MDSplus data.

beta_N = NaN(size(timebase));
beta_p = NaN(size(timebase));
dbeta_p_dt = NaN(size(timebase));
kappa = NaN(size(timebase));
upper_gap = NaN(size(timebase));
lower_gap = NaN(size(timebase));
li = NaN(size(timebase));
dli_dt = NaN(size(timebase));
q0 = NaN(size(timebase));
qstar = NaN(size(timebase));
q95 = NaN(size(timebase));
V_loop_efit = NaN(size(timebase));
Wdiam = NaN(size(timebase));
dWdiam_dt = NaN(size(timebase));

% All the signals used inside this routine are column vectors, since
% mdsvalue returns column vectors.  But the input parameter, 'timebase',
% may be either a column or row vector.  Therefore, in order to continue
% with our column-based calculations, we need to create a copy of
% 'timebase' that is guaranteed to be a column vector.  Just before
% returning from this routine, all output vectors will be converted to the
% same shape (i.e. column or row vector) as the 'timebase' input
% parameter.

if (size(timebase,1) > 1);
  timebase_column = timebase;
else;
  timebase_column = transpose(timebase);
end;

% EFIT runs that are done automatically after each shot are usually put in
% the EFIT01 and EFIT02 trees.  Additional EFIT runs using non-standard
% settings and/or time sampling are often done by individuals, and these go
% into trees named EFIT03, EFIT04, etcetera.  This includes specialized
% EFIT runs that we have done for use with our disruption warning database.
% An SQL database of available EFIT trees exists, from which our specific
% EFIT tree can be selected by means of the 'run_by' and/or 'runtype'
% keywords.  All of our EFIT runs have runtype = 'DIS' (for disruption).

efittrees  = select_efit_trees(shot, '', 'DIS');
%efittrees  = select_efit_trees(shot, 'granetzr', 'DIS');
%efittrees = select_efit_trees(shot, {'granetzr', 'reac'}, 'DIS');

if isempty(efittrees);
  fprintf(1, 'No disruption EFIT tree for this shot\n');
  return;
end;
tree=char(efittrees(end,:));
[shotopened, status] = mdsopen(tree, shot);
if (mod(status,2) == 0);
  return;
end;

% Read in EFIT timebase.  DIII-D is "stupid" (Bob quote) and records time
% in ms. Also, note the efit timebase data is in a node called "atime"
% instead of "time" (where "time" does not work).

[efittime, status] = mdsvalue('\efit_a_eqdsk:atime');
efittime = efittime/1000; % efit time in seconds

if (mod(status,2)==0 || length(efittime) <= 4);
  return;
end;

beta_N = mdsvalue('\efit_a_eqdsk:betan');
beta_p = mdsvalue('\efit_a_eqdsk:betap');
kappa = mdsvalue('\efit_a_eqdsk:kappa'); % changed from "eout" from C-Mod
li = mdsvalue('\efit_a_eqdsk:li');
upper_gap = mdsvalue('\efit_a_eqdsk:gaptop'); % meters, changed from "otop"
                                             % from C-Mod 
lower_gap = mdsvalue('\efit_a_eqdsk:gapbot'); % meters, changed from
                                              % "obott" from C-Mod 
q0 = mdsvalue('\efit_a_eqdsk:q0');
qstar = mdsvalue('\efit_a_eqdsk:qstar');
q95 = mdsvalue('\efit_a_eqdsk:q95');
V_loop_efit = mdsvalue('\efit_a_eqdsk:vsurf'); % Carlos Paz-Soldan
                                               % recommended that we use the
					       % vsurf signal
Wdiam = mdsvalue('\efit_a_eqdsk:wmhd'); % wmhd is the EFIT reconstruction
                                        % of stored energy 
chisq = mdsvalue('\efit_a_eqdsk:chisq'); % Use chisq to determine which time
                                         % slices are invalid.
mdsclose;

% EFIT reconstructions are sometimes invalid, particularly when very close
% to a disruption.  There are a number of EFIT parameters that can indicate
% invalid reconstructions, such as 'terror' and 'chisq'.  Here we use
% 'chisq' to determine which time slices should be excluded from our
% disruption warning database.

invalid_indices = find(chisq > 50);

beta_N(invalid_indices) = NaN;
beta_p(invalid_indices) = NaN;
kappa(invalid_indices) = NaN;
li(invalid_indices) = NaN;
upper_gap(invalid_indices) = NaN;
lower_gap(invalid_indices) = NaN;
q0(invalid_indices) = NaN;
qstar(invalid_indices) = NaN;
q95(invalid_indices) = NaN;
V_loop_efit(invalid_indices) = NaN;
Wdiam(invalid_indices) = NaN;

dbeta_p_dt = gradient(beta_p, efittime);
dli_dt = gradient(li, efittime);
dWdiam_dt = gradient(Wdiam, efittime);

% Interpolate data onto the requested timebase

beta_N = interp1(efittime, beta_N, timebase_column, 'linear');
beta_p = interp1(efittime, beta_p, timebase_column, 'linear');
kappa = interp1(efittime, kappa, timebase_column, 'linear');
li = interp1(efittime, li, timebase_column, 'linear');
upper_gap = interp1(efittime, upper_gap, timebase_column, 'linear');
lower_gap = interp1(efittime, lower_gap, timebase_column, 'linear');
q0= interp1(efittime, q0, timebase_column, 'linear');
qstar = interp1(efittime, qstar, timebase_column, 'linear');
q95 = interp1(efittime, q95, timebase_column, 'linear');
V_loop_efit = interp1(efittime, V_loop_efit, timebase_column, 'linear');
Wdiam = interp1(efittime, Wdiam, timebase_column, 'linear');
dbeta_p_dt = interp1(efittime, dbeta_p_dt, timebase_column, 'linear');
dli_dt = interp1(efittime, dli_dt, timebase_column, 'linear');
dWdiam_dt = interp1(efittime, dWdiam_dt, timebase_column, 'linear');

if (size(timebase,2) > 1);    % If "timebase" is a row vector, then convert
  beta_N = transpose(beta_N); % all the outputs back to row vectors
  beta_p = transpose(beta_p);
  kappa = transpose(kappa);
  li = transpose(li);
  upper_gap = transpose(upper_gap);
  lower_gap = transpose(lower_gap);
  q0 = transpose(q0);
  qstar = transpose(qstar);
  q95 = transpose(q95);
  V_loop_efit = transpose(V_loop_efit);
  Wdiam = transpose(Wdiam);
  dbeta_p_dt = transpose(dbeta_p_dt);
  dli_dt = transpose(dli_dt);
  dWdiam_dt = transpose(dWdiam_dt);
end;

end
