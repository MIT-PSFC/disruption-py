function prad_peaking = get_prad_peaking(shot, timebase);

% This routine calculates the peaking factor of the profiles of radiated
% power measured by the AXUV arrays on EAST.  Here we define the peaking
% factor as the ratio of the average of the Prad signals of the core
% channels to the average of all the channels, excluding those that look in
% the divertor region (core-to-average).  This routine linearly
% interpolates the prad_peaking signal onto the given timebase.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   prad_peaking = ratio of core Prad signals to all Prad signals
%                  Note: for now we are defining "core" to be the
%                  centralmost 6 chords out of the 64 chords in the EAST
%                  axuv arrays, and "all" to be all of the
%                  non-divertor-viewing chords.  See comments later in this
%                  code for more details.
%
% Author: Robert Granetz   2019/05/13
%
% Revisions:
%  2019/mm/dd -- R. Granetz: 

% The following section about calibration factors was copied-and-pasted
% directly from Duan Yanmin <ymduan@ipp.ac.cn>'s code.

%%%%%%%%%%%%calibration factors%%%
Fac1=[1.3681,1.3429,1.3215,1.3039,1.2898,1.2793,1.2723,1.2689,1.2689,1.2723,1.2793,1.2898,1.3039,1.3215,1.3429,1.3681]*10.^(4);
Fac2=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
      1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];% factors of Amp.Gain
Fac3=[1,1,1,1]; %cross calibration factors between arrays
Fac4=1*10.^(-3);%unit convert
Fac5=2.5;% corrected factor by cross calibration with foil bolometer
Maj_R=1.85;
Del_r=[3.6,3.6,3.5,3.4,3.4,3.3,3.3,3.2,3.2,3.1,3.1,3.0,3.0,2.9,2.9,2.8,...
  2.9,2.8,2.8,2.8,2.8,2.8,2.8,2.7,2.7,2.7,2.7,2.7,2.6,2.6,2.6,2.6,...
  2.6,2.6,2.6,2.6,2.7,2.7,2.7,2.7,2.7,2.8,2.8,2.8,2.8,2.8,2.8,2.9,...
  2.8,2.9,2.9,3.0,3.0,3.1,3.1,3.2,3.2,3.3,3.3,3.4,3.4,3.5,3.6,3.6]*0.01;

%%%%%%%%%%%%%%%%%

prad_peaking = NaN(size(timebase));

mdsconnect('mds.ipp.ac.cn');
[~, status] = mdsopen('east_1', shot);
if (mod(status, 2) == 0);
 return;
end;

[xuvtime, status] = mdsvalue('dim_of(\pxuv1)');
if (mod(status, 2) == 0);
  mdsclose;
  return;
end;

ntimes = length(xuvtime);
xuv = NaN(ntimes, 64);  % There are 64 AXUV chords, arranged in 4 arrays
                        % of 16 channels each.

for iarray = 1:4;
  for ichan = 1:16;
    ichord = (iarray-1) * 16 + ichan;
    [signal, status] = mdsvalue(['\pxuv' num2str(ichord, '%i')]);
    xuv(:, ichord) = signal;
    xuv(:, ichord) = xuv(:,ichord) - mean(xuv(1:100, ichord)); % substract
                                                               % baseline
    dt = xuvtime(2) - xuvtime(1); % sampling interval in seconds
    smoothing_time = 1.e-3;
    smoothing_window = max([round(smoothing_time/dt), 1]);
    xuv(:, ichord) = smooth(xuv(:, ichord), smoothing_window);
    xuv(:, ichord) = xuv(:, ichord) * ...
      Fac1(ichan) * Fac2(iarray, ichan) * Fac3(iarray) * Fac4 * ...
      2*pi * Maj_R * Del_r(ichan);  % from Duan Yanmin's program
  end
end
mdsclose;

% correction for bad channels (from Duan Yanmin's program)
xuv(:, 12) = 0.5 * (xuv(:, 11) + xuv(:, 13));
%xuv(:, 36) = 0.5 * (xuv(:, 35) + xuv(:, 36));

% Define the core chords to be #28 to #37 (centermost 10 chords), and
% define the non-divertor chords to be #09 to #56.  (Chords #1-8 view the
% lower divertor region, and chords #57-64 view the upper divertor region.)

core = zeros(64,1);
core(30:35) = 1/length([30:35]);
all = zeros(64,1);
all(9:56) = 1/length([9:56]);

prad_peaking = NaN(size(xuvtime));

for itime = 1:ntimes;
  prad_peaking(itime) = (xuv(itime, :) * core) / (xuv(itime, :) * all);
end;
prad_peaking = interp1(xuvtime, prad_peaking, timebase);

end
