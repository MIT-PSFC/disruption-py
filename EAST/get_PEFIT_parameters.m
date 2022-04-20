function [beta_n, beta_p, kappa, li, q95, Wmhd, aminor] = ...
  get_PEFIT_parameters(shot, timebase);

% This function obtains the P-EFIT parameters from 'pefit_east' tree
%
% Inputs:
%   shot = shot number
%   timebase = array of times at which to evaluate the parameters
%
% Outputs:
%   beta_n = normalized beta
%   beta_p = beta poloidal
%   kappa = elongation
%   li = internal inductance
%   q95 = edge safety factor
%   Wmhd = diamagnetic energy, [J]
%   aminor = plasma minor radius, [m]
%  
%
% Author: Cristina Rea 2019/05/13 get_EFIT_parameters_EAST was
%	  adapted to retrieve only a reduced number of parameters
%	  from P-EFIT runs (Yao Huang reran these for us).
%
% Revision/history:
%  2017/12/13 - RSG; Use chisq to determine invalid time slices

% Initialize all outputs to NaN (Not-a-Number).  These will be returned by
% this routine if it is unable to open the MDSplus tree or read out the
% MDSplus data.

beta_n     = NaN(size(timebase));
beta_p     = NaN(size(timebase));
kappa      = NaN(size(timebase));
li         = NaN(size(timebase));
q95        = NaN(size(timebase));
Wmhd       = NaN(size(timebase));
aminor     = NaN(size(timebase));

% The input array, "timebase", can be either a row vector (1 x n) or a
% column vector (n x 1).  We desire to have the output arrays match the
% shape of "timebase".  However, in Matlab the "mdsvalue" routine returns
% column vectors for 1-D signals, so it is simpler to work with column
% vectors within this routine, and then, if "timebase" is a row vector,
% convert the outputs to row vectors just before exiting this routine.  So
% the first step is to create a copy of "timebase" that is guaranteed to be
% a column vector.

if (size(timebase,1) > 1);
  timebase_column = timebase;
else;
  timebase_column = transpose(timebase);
end;

% Open 'PEFIT_EAST' tree.  If error, try opening the standard EFIT tree,
% 'EFIT_EAST'.  If that also fails, then exit.

[shotopened, status] = mdsopen('pefit_east', double(shot));
if (mod(status,2) == 0);
  [shotopened, status] = mdsopen('efit_east', double(shot));
  if (mod(status,2) == 0);
    return;
  end;
end;
% Read in the EFIT timebase.  If error, exit routine.

[efittime, status] = mdsvalue('\efit_aeqdsk:atime');
[efittime, unique_indices, ~] = unique(efittime);
if (mod(status,2)==0 || length(efittime) <= 1);
  return;
end;

% Read in the P-EFIT signals.  Assume that if the first signal is okay, all
% the signals will all be okay.  And if the first signal has an error, they
% will all have errors, and we should exit the routine.

[beta_n, status] = mdsvalue('\efit_aeqdsk:betan');
if (mod(status,2)==0);
  beta_n = NaN(size(timebase));
  return;
end;
beta_p = mdsvalue('\efit_aeqdsk:betap'); 
kappa = mdsvalue('\efit_aeqdsk:kappa'); % different from C-Mod
li = mdsvalue('\efit_aeqdsk:li');
q95 = mdsvalue('\efit_aeqdsk:q95');
Wmhd = mdsvalue('\efit_aeqdsk:wmhd');
aminor = mdsvalue('\efit_aeqdsk:aminor');
chisq = mdsvalue('\efit_aeqdsk:chisq'); % Use chisq to determine which time
                                        % slices are invalid.
convergence = mdsvalue('\efit_aeqdsk:error'); % Use convergence error to determine
                                              % invalid time slices.
mdsclose;
  
beta_n = beta_n(unique_indices);
beta_p = beta_p(unique_indices);
kappa = kappa(unique_indices);
li = li(unique_indices);
q95 = q95(unique_indices);
Wmhd = Wmhd(unique_indices);
aminor = aminor(unique_indices);
chisq = chisq(unique_indices);
% convergence = convergence(unique_indices);

% P-EFIT reconstructions are sometimes invalid, particularly when very close
% to a disruption.  There are a number of P-EFIT parameters that can indicate
% invalid reconstructions, such as 'error' and 'chisq'.  Here we use
% 'chisq' to determine which time slices should be excluded from our
% disruption warning database.
% Yao Huang suggests to use data with:
%   - chisq < 20
%   - convergence error < 1 (data stored in MDS+ is multiplied by 1e3)
%   - ip > 180 kA
% For now, we only check the first two conditions. 
% If ever we want to extend analysis to ramp up or down we need to check ip.

invalid_indices = find(chisq > 20 & convergence < 1);

beta_n(invalid_indices) = NaN;
beta_p(invalid_indices) = NaN;
kappa(invalid_indices) = NaN;
li(invalid_indices) = NaN;
q95(invalid_indices) = NaN;
Wmhd(invalid_indices) = NaN;
aminor(invalid_indices) = NaN;

% Interpolate data onto the requested timebase.  This routine used to have
% 'interp1' extrapolate, using the value 'zero' for the extrapolation.  We
% should not do that for EFIT parameters.  If extrapolation is needed, the
% extrapolated values should be 'NaN' for EFIT parameters.

beta_n    = interp1(efittime, beta_n,    timebase_column);
beta_p    = interp1(efittime, beta_p,    timebase_column);
kappa     = interp1(efittime, kappa,     timebase_column);
li        = interp1(efittime, li,        timebase_column);
q95       = interp1(efittime, q95,       timebase_column);
Wmhd      = interp1(efittime, Wmhd,      timebase_column);
aminor    = interp1(efittime, aminor,    timebase_column);

% If "timebase" is a row vector, then convert all the outputs back to row
% vectors before exiting this routine.

if (size(timebase,2) > 1);
  beta_n    = transpose(beta_n);
  beta_p    = transpose(beta_p);
  kappa     = transpose(kappa);
  li        = transpose(li);
  q95       = transpose(q95);
  Wmhd      = transpose(Wmhd);
  aminor    = transpose(aminor);
end;

end
