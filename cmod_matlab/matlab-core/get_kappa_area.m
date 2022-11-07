function kappa_area = get_kappa_area(shot, timebase);

% This function obtains the EFIT parameters, and also calculates the time
% derivatives of a few of them.
%
% Inputs:
%   shot = shot number
%   timebase = array of times at which to evaluate the parameters
%
% Outputs:
%   kappa_area = elongation
%
% Author: C. Rea 2020-03-10 adapted from get_EFIT_parameters to get kappa_area.
%                           This is an alternate definition of elongation:
%                            kappa_area = (plasma cross-section area)/pi*a^2
%                           Note: 'cross-section area' and 'a' come from EFIT.
% 2020/03/11 - RSG; modified to work on Alcator C-Mod EFIT data

% Initialize all outputs to NaN (Not-a-Number).  These will be returned by
% this routine if it is unable to open the MDSplus tree or read out the
% MDSplus data.

kappa_area = NaN(size(timebase));

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

% Read in EFIT timebase

[efittime, status] = mdsvalue('\efit_aeqdsk:time');
if (mod(status,2)==0 || length(efittime) == 1);
  return;
end;

aminor = mdsvalue('\efit_a_eqdsk:aminor');
area = mdsvalue('\efit_a_eqdsk:area');

mdsclose;

% compute kappa_area
kappa_area = area./(pi*aminor.^2);

% Interpolate data onto the requested timebase
kappa_area = interp1(efittime, kappa_area, timebase_column, 'linear');

if (size(timebase,2) > 1);    % If "timebase" is a row vector, then convert
  kappa_area = transpose(kappa_area); % all the outputs back to row vectors
end;

end
