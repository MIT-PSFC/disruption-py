function [kappa_area] = get_kappa_area(shot, timebase);

% This function obtains kappa_area from EFIT parameters

% Inputs:
%   shot = shot number
%   timebase = array of times at which to evaluate the parameters
%
% Outputs:
%   kappa_area = elongation (from cross section area)
%  
% Author: Modified on EAST by Wang Bo 2015/12/10
%         Alex Tinguely 2015-10-07 
%         Some modifications on EAST by Robert Granetz 2016/04/05
%         2017/04/10 - R. Granetz; modified so that program tries to the
%           standard efit tree, 'efit_east', in the case when our special
%           efit tree, 'efit18', does not exist.
%         Modified to have fewer data by C. Rea 2018/08/08
%         Mod to get kappa_area by C. Rea 2020/03/11

% Revision/history:
%  2017/12/13 - RSG; Use chisq to determine invalid time slices

% Initialize all outputs to NaN (Not-a-Number).  These will be returned by
% this routine if it is unable to open the MDSplus tree or read out the
% MDSplus data.

kappa_area       = NaN(size(timebase));

fprintf(1,['Processing kappa_area for shot %7i\n'], shot);
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

% Open EFIT18 tree.  If error, try opening the standard EFIT tree,
% 'EFIT_EAST'.  If that also fails, then exit.

[shotopened, status] = mdsopen('efit18', double(shot));
if (mod(status,2) == 0);
  [shotopened, status] = mdsopen('efit_east', double(shot));
  if (mod(status,2) == 0);
    fprintf(1,'no efit_east \n');
    return;
  end;
end;
% Read in the EFIT timebase.  If error, exit routine.

[efittime, status] = mdsvalue('\efit_aeqdsk:atime');
[efittime, unique_indices, ~] = unique(efittime);
if (mod(status,2)==0 || length(efittime) <= 1);
  return;
end;

% Read in the EFIT signals.  Assume that if the first signal is okay, all
% the signals will all be okay.  And if the first signal has an error, they
% will all have errors, and we should exit the routine.

% I need the minor radius, aminor, to calculate kappa_area. 
% I get this from EFIT (\efit_aeqdsk:aout)

% Open EFIT18 tree.  If error, try opening the standard EFIT tree,
% 'EFIT_EAST'.  If that also fails, then use a default value of 0.45 m

[aminor, status] = mdsvalue('\efit_aeqdsk:aout');
if (mod(status,2) == 0);
  aminor = 0.45 * ones(size(efittime));
end;
[area, status] = mdsvalue('\efit_aeqdsk:area');
if (mod(status,2)==0);
  kappa_area = NaN(size(timebase));
  return;
end;

chisq = mdsvalue('\efit_aeqdsk:chisq'); % Use chisq to determine which time
                                        % slices are invalid.
mdsclose;
  
area = area(unique_indices);
aminor = aminor(unique_indices);
chisq = chisq(unique_indices);

% calculate kappa_area, convert area from cm^2 to m^2
kappa_area = area*1e-4 ./ (pi*aminor.^2);

% EFIT reconstructions are sometimes invalid, particularly when very close
% to a disruption.  There are a number of EFIT parameters that can indicate
% invalid reconstructions, such as 'terror' and 'chisq'.  Here we use
% 'chisq' to determine which time slices should be excluded from our
% disruption warning database.

invalid_indices = find(chisq > 20);

kappa_area(invalid_indices) = NaN;

% Interpolate data onto the requested timebase.  This routine used to have
% 'interp1' extrapolate, using the value 'zero' for the extrapolation.  We
% should not do that for EFIT parameters.  If extrapolation is needed, the
% extrapolated values should be 'NaN' for EFIT parameters.

kappa_area     = interp1(efittime, kappa_area,      timebase_column);

% If "timebase" is a row vector, then convert all the outputs back to row
% vectors before exiting this routine.

if (size(timebase,2) > 1);
  kappa_area     = transpose(kappa_area);
end;

end
