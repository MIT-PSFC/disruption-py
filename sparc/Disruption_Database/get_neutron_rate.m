function rate = get_neutron_rate(shot, timebase)

% This script gets the neutron rate in neutrons/s.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   rate = neutron rate (1/s)
%
% Author: Alex Tinguely 2015-09-09

[timebase, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

mdsopen('particles', shot);

Rate = mdsvalue('\PARTICLES::TOP.NEUTRONS.HE_3_BANK.RESULTS:HE3_NRATE'); % actual/measured neutrons/s
t = mdsvalue('dim_of(\PARTICLES::TOP.NEUTRONS.HE_3_BANK.RESULTS:HE3_NRATE)'); % time

if strcmp(Rate(1), 'J') || strcmp(t(1), 'J') || length(Rate) == 1 || length(t) == 1
    % if Rate,t are not callable or only one point
    rate = zeros(1, length(timebase)) + NaN; % return NaN
    return 
end

rate = interp1(t, Rate, timebase); % interpolate the neutron rate over given timebase

mdsclose;

% reorient all output arrays to original configuration
rate = orient_array(rate, original_orientation);

end