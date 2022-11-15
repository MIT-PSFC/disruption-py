function [kappa_area] = get_kappa_area(shot, timebase);

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
% Author: Alex Tinguely 2015-10-07 (Dan Brunner's 30th Birthday!)
% Some modifications by Robert Granetz 2015/10/20
% Updates: DIII-D by Alex Tinguely 2016-02-26
% Updates: DIII-D EFIT parameters changes by Cristina Rea 2016-08-11
% Updates: C. Rea 2020-03-10 modified to get kappa_area

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

% EFIT runs that are done automatically after each shot are usually put in
% the EFIT01 and EFIT02 trees.  Additional EFIT runs using non-standard
% settings and/or time sampling are often done by individuals, and these go
% into trees named EFIT03, EFIT04, etcetera.  This includes specialized
% EFIT runs that we have done for use with our disruption warning database.
% An SQL database of available EFIT trees exists, from which our specific
% EFIT tree can be selected by means of the 'run_by' and/or 'runtype'
% keywords.  All of our EFIT runs have runtype = 'DIS' (for disruption).

efittrees  = select_efit_trees(shot, '', 'DIS');

if isempty(efittrees);
  fprintf(1, 'No disruption EFIT tree for this shot\n');
  fprintf(1, 'Using EFIT01\n');
  efittrees = {'EFIT01'};
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
efittime = efittime*1e-3; % efit time in seconds

if (mod(status,2)==0 || length(efittime) <= 4);
  return;
end;

aminor = mdsvalue('\efit_a_eqdsk:aminor');
area = mdsvalue('\efit_a_eqdsk:area');
chisq = mdsvalue('\efit_a_eqdsk:chisq'); % Use chisq to determine which time
                                         % slices are invalid.
mdsclose;

% compute kappa_area
kappa_area = area./(pi*aminor.^2);

% EFIT reconstructions are sometimes invalid, particularly when very close
% to a disruption.  There are a number of EFIT parameters that can indicate
% invalid reconstructions, such as 'terror' and 'chisq'.  Here we use
% 'chisq' to determine which time slices should be excluded from our
% disruption warning database.

invalid_indices = find(chisq > 50);

kappa_area(invalid_indices) = NaN;

% Interpolate data onto the requested timebase
kappa_area = interp1(efittime, kappa_area, timebase_column, 'linear');

if (size(timebase,2) > 1);    % If "timebase" is a row vector, then convert
  kappa_area = transpose(kappa_area); % all the outputs back to row vectors
end;

end
