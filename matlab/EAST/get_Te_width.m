function [Te_width_normalized,Te_width_MI] = get_Te_width(shot, timebase);

% This routine calculates a representative half-width of the Te profile
% measured by the Michelson interferometer on EAST.  The profiles extend
% from the center to the outboard edge of the plasma.  The diagnostic
% measurement is corrupted by non-thermal emission due to lower hybrid,
% especially for R > 2.15 m.  Unfortunately most shots on EAST have LH
% power turned on for most of the discharge.  So I only use the portion of
% the profiles from R = 1.8 to 2.15 m.  Over this restricted radial range,
% the profiles typically have noticeable wiggles, even though the signals
% have all been calibrated.  But for now, this is the best we have.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Output:
%   Te_width_MI = Te profile width from Michelson interferometer [m]
%   Te_width_normalized = Te_width_MI normalized by plasma minor radius a
%
% Author: Robert Granetz   June 2016
%
% Note: the person responsible for the Michelson interferometer
%       is Liu Yong <liuwong@ipp.ac.cn>

% The input array, "timebase", can be either a row vector (1 x n) or a
% column vector (n x 1).  We desire to have the output arrays match the
% shape of "timebase".  However, the "mdsvalue" routine in Matlab returns
% column vectors for 1-D signals, so it is simpler to work with column
% vectors within this routine, and then, if "timebase" is a row vector,
% convert the outputs to row vectors just before exiting this routine.  So
% the first step is to create a copy of "timebase" that is guaranteed to be
% a column vector.
%
% Revision/history:
%  2018/10   - Chen Dalong, K. Montes, RSG: added Te_width_normalized
%              and associated code

if (size(timebase,1) > 1);
  timebase_column = timebase;
else;
  timebase_column = transpose(timebase);
end;

% Initialize all output arrays to NaN (Not-a-Number) column vectors

Te_width_MI = NaN(length(timebase), 1);

% Read in the Michelson interferometer data from MDSplus

[shotopened, status] = mdsopen('analysis', double(shot));
if (mod(status,2) == 1);
  [Te, status] = mdsvalue('\te_mi');  % Read in Te(r,t) [eV]
  if (mod(status,2) == 1);            % If successful, continue
    R = mdsvalue('\r_mi');
    time = mdsvalue('dim_of(\te_mi,0)');

% Use only data in the range 1.80 <= R < 2.15 m

    rindx = find(R >= 1.80 & R < 2.15);
    R = R(rindx);
    Te = Te(rindx,:);

% For each time slice, fit the Te profiles with a parabola (polynomial fit
% of degree = 2), and interpolate the fit onto a dense radial coordinate,
% and find the half-width at half max (HWHM).

    rarray = 1.45 : .001 : 2.15;
    Te_width_MI = NaN(length(time), 1);
    for it = 1:length(time);
      if (length(R) < 4 || length(Te) < 4);
        continue;
      end;
      p = polyfit(R, Te(:, it), 2);
      Te_array = polyval(p, rarray);
      [Te_max, maxindx] = max(Te_array);
      r_max = rarray(maxindx);
      Te_HM = Te_max/2;
      [~, HM_indices] = min(abs(Te_array - Te_HM));
      r_HM = rarray(HM_indices(1));
      Te_width_MI(it) = abs(r_max - r_HM);
    end;
    Te_width_MI = interp1(time, Te_width_MI, timebase_column);
  end;
  mdsclose;
end;

% Retreve plasma minor radius 'aminor' from the EFIT tree and use it to
% normalize the temperature profile width

[shotopened, status] = mdsopen('efit_east', double(shot));
if (mod(status,2) == 1);
  [aminor, status] = mdsvalue('\aminor');  % Read in Te(r,t) [eV]
  if (mod(status,2) == 1);                 % If successful, continue
    aminortime = mdsvalue('dim_of(\aminor)');
    if (length(aminortime) >= 2);  % deal with rare bug
      aminor = interp1(aminortime,aminor,timebase_column);
    else;
      aminor = NaN(size(Te_width_MI));
    end;
  else
    aminor = NaN(size(Te_width_MI));
  end
  mdsclose;
else
  aminor = NaN(size(Te_width_MI));
end

Te_width_normalized = Te_width_MI./aminor;

% The output signal(s) are currently column vectors.  However, we desire to
% have the output array(s) match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase,2) > 1);
  Te_width_MI = transpose(Te_width_MI);
  Te_width_normalized = transpose(Te_width_normalized);
end;

end
