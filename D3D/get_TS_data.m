function [Te_HWHM, Te_width_normalized] = get_TS_data( shot, timebase);

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

if (size(timebase,1) > 1);
  timebase_column = timebase;
else;
  timebase_column = transpose(timebase);
end;

mdsconnect('atlas.gat.com');

[shotopened, status]=mdsopen('electrons', shot);
if (mod(status,2)==0);
% fprintf(1,'  Unable to open ELECTRONS tree for shot%7i\n', shot);
  Te_HWHM = NaN(size(timebase));
  Te_width_normalized = NaN(size(timebase));
  return;
end;

% Read in Thomson core temperature data, which is a 2-D array, with the
% dependent dimensions being time and z (vertical coordinate)

[TS_time, status] = mdsvalue('dim_of(\top.ts.blessed.core:temp,0)');

if (mod(status,2) == 1 & length(TS_time) >= 2);
  TS_data = mdsvalue('\top.ts.blessed.core:temp');
  TS_z = mdsvalue('dim_of(\top.ts.blessed.core:temp,1)');
  TS_time = TS_time/1.e3; % convert ms to s
% Get rid of the last channel (#41), which is not real
  TS_data = TS_data(:, 1:end-1);
  TS_z = TS_z(1:end-1);
  mdsclose;

  Te_HWHM = NaN(length(TS_time),1);

  zarray = [0:.01:.9];
  itimes = find(TS_time > 0);

  for i = 1:length(itimes);
    y = TS_data(itimes(i),:);
    ok_indices = find(y ~= 0);
    y = y(ok_indices);
    z = TS_z(ok_indices);
    if (length(ok_indices) > 2);
      p = polyfit(z, transpose(y), 2);
      Te_array = polyval(p, zarray);
      [Te_max, maxindx] = max(Te_array);
      z_max = zarray(maxindx);
      Te_HM = Te_max/2;
      [~, HM_indices] = min(abs(Te_array - Te_HM));
      HM_indx = max(HM_indices);
      z_HM = zarray(HM_indx);
      if z_HM > z_max;
        Te_HWHM(itimes(i)) = z_HM - z_max;
      end;
    end;
  end;
  Te_HWHM = interp1(TS_time, Te_HWHM, timebase_column, 'linear');
else;
  mdsclose;
  Te_HWHM = NaN(size(timebase));
  Te_width_normalized = NaN(size(timebase));
  return;
end;

% Get minor radius from EFIT.  This is needed to normalize the Te width.

efittrees = select_efit_trees(shot,'granetzr','DIS');
if isempty(efittrees);
% fprintf(1, 'No disruption EFIT tree for this shot\n');
  [shotopened, status] = mdsopen('efit01', shot);
  if (mod(status,2) == 0);
%   return;
  end;
else;
  tree=char(efittrees(end,:));
  [shotopened, status] = mdsopen(tree, shot);
  if (mod(status,2) == 0);
%   return;
  end;
end;

% Read in EFIT minor radius and timebase.  However, if the minor radius
% data is not available, use a default fixed value of 0.59 m.  (We surveyed
% several hundred shots to determine this default value.)  Also, note the
% efit timebase data is in a node called "atime" instead of "time" (where
% "time" does not work).

[efittime, status] = mdsvalue('\efit_a_eqdsk:atime');

if (mod(status,2)==1 && length(efittime) > 4);
  efittime = efittime/1.e3; % convert ms to s
  aminor = mdsvalue('\efit_a_eqdsk:aminor'); % meters
  mdsclose;
% Interpolate data onto the requested timebase
  aminor = interp1(efittime, aminor, timebase_column, 'linear');
  Te_width_normalized = Te_HWHM ./ aminor;
else;
  mdsclose;
  Te_width_normalized = Te_HWHM / 0.59;
end;

% The output signals are currently all column vectors.  However, we desire
% to have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase,2) > 1);
  Te_HWHM = transpose(Te_HWHM);
  Te_width_normalized = transpose(Te_width_normalized);
end;

end
