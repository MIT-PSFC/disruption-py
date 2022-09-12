function peaking_factor = get_pressure_peaking_factor(shot, timebase)

% This function calculates the pressure peaking factor which is the ratio
% of the on-axis pressure to the average pressure. The on-axis pressure p0
% is calculated using the electron temperature and density on-axis. The
% average pressure is calculated using the plasma energy and volume stored
% in EFIT.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   peaking_factor = p0/p_avg = p0/<p>
%
% Author: Alex Tinguely 15-09-16

[timebase, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

mdsopen('electrons', shot);

e = 1.6021766208e-19; % 1eV = e Joules, also charge of electron

T = mdsvalue('\electrons::thom_midpln:te_t'); % [keV], electron temperature on axis in keV
Te0 = T*1e3*e; % [J], electron temperature on axis in Joules

t0 = mdsvalue('dim_of(\electrons::thom_midpln:te_t)'); % time, should be same for ne0
ne0 = mdsvalue('\electrons::thom_midpln:ne_t'); % [1/m^3], electron density on axis

if strcmp(T(1), 'J') || strcmp(ne0(1), 'J') % if either Te0 or ne0 is not callable (there is an error)
    % so the first letter of the error is 'J'
    p0 = zeros(1,length(timebase)) + NaN; % set p0 to NaN
    
else   
    
    p0 = interp1(t0, 2.*ne0.*Te0, timebase); % [J/m^3 = N/m^2 = Pa] 
    % p0 = 2*ne0*Te0 is upper bound on pressure
    % becasue we are assuming that ne0 = ni0 and Te0 = Ti0 (i = ions) when we usually
    % have ne0 >= ni0 and Te0 >= Ti0
    
end

mdsclose;

mdsopen('mhd', shot);

vol = mdsvalue('\ANALYSIS::EFIT_AEQDSK:VOLUME'); % plasma volume (m^3)
t = mdsvalue('dim_of(\ANALYSIS::EFIT_AEQDSK:VOLUME)'); % EFIT time, should be same for Wp
Wp = mdsvalue('\ANALYSIS::EFIT_AEQDSK:WPLASM'); % plasma stored energy [J]

p_avg = interp1(t, 2/3.*Wp./vol, timebase); % average pressure = <p> = 2/3*energy/volume 
                                            % from energy = 3/2*<pressure>*volume (ideal gas laws)

peaking_factor = p0./p_avg; % ratio of on-axis to average pressure, = p0/<p>

mdsclose;

% reorient all output arrays to original configuration
peaking_factor = orient_array(peaking_factor, original_orientation);

end