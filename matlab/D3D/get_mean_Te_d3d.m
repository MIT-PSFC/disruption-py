function Te = get_mean_Te_d3d( shot, timebase)

% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Output:
%   Te_HWHM(timebase) = half-width at half-max of parabolic fits to the
%     Te(t,z) data
%   Te_width_normalized(timebase) = half-width normalized to the minor
%     radius, which is obtained from EFIT
%
% Authors: Robert Granetz and Cristina Rea        May 2017
% Revision history:
%  2018/01/29 -- R. Granetz; added Te_width_normalized

% The input array, "timebase", can be either a row vector (1 x n) or a
% column vector (n x 1).  We desire to have the output arrays match the
% shape of "timebase".  However, Matlab calls to the routine "mdsvalue"
% return column vectors for 1-D signals, so it is simpler to work with
% column vectors within this routine, and then, if necessary, convert the
% outputs to row vectors just before exiting this routine.  So the first
% step is to create a copy of "timebase" that is guaranteed to be a column
% vector.

if (size(timebase,2) > 1)
  timebase = transpose(timebase);
end
Te = NaN(size(timebase));
mdsconnect('atlas.gat.com');

[~, status]=mdsopen('electrons', shot);
if (mod(status,2)==0)
% fprintf(1,'  Unable to open ELECTRONS tree for shot%7i\n', shot);
  mdsclose;
  return
end

% Read in Thomson core temperature data, which is a 2-D array, with the
% dependent dimensions being time and z (vertical coordinate)

[TS_time, status] = mdsvalue('dim_of(\top.ts.blessed.core:temp,0)');

if (mod(status,2) == 1 && length(TS_time) >= 2)
  TS_data = mdsvalue('\top.ts.blessed.core:temp');
  
  TS_time = TS_time/1.e3; % convert ms to s
% Get rid of the last channel (#41), which is not real
  TS_data = TS_data(:, 1:end-1);
  ind = find(TS_data<0);
  TS_data(ind) = NaN;
  mdsclose;
  indx= TS_time>0;
  TS_time=TS_time(indx);
  if isempty(TS_time)
      return
  end
  

  Te_0 = mean(TS_data,2,'omitnan');
  Te_0 = Te_0(indx);
  id = find(~isnan(Te_0), 1); 
  if isempty(id) || isempty(Te_0)
        return
  end
  Te = interp1(TS_time,Te_0,timebase);
else
    mdsclose;
    return
end

  