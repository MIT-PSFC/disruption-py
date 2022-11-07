function [H89] = get_H89(shot, timebase)

% This function calculates H89, the ratio between the energy confinement
% time, tE, and ITER-89 L-Mode Scaling Law for confinement time, t89.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   H89 = tE/t89
%
% Author: Alex Tinguely 2015-09-18

R0 = 0.68; %[m] C-Mod major radius from wiki
a = 0.22; % [m] C-Mod minor radius from wiki

[timebase, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

H89 = zeros(1, length(timebase)) + NaN; % initialize to NaN

mdsopen('mhd', shot);
    
    wp = mdsvalue('\ANALYSIS::EFIT_AEQDSK:WPLASM'); % plasma stored energy [J]
    t = mdsvalue('dim_of(\ANALYSIS::EFIT_AEQDSK:WPLASM)'); % EFIT time
    
    if strcmp(wp(1), 'J') || strcmp(t(1), 'J') || length(wp) == 1 || length(t) == 1
        % if wp,t are not callable or only one point
        return 
    end
    
    Wp = interp1(t, wp, timebase); % interpolated energy [J]

    % if the plasma energy is less than zero (sometimes it is?), set it to zero
    for i = 1:length(timebase)
       if Wp(i) < 0
          Wp(i) = 0;
       end
    end
    
    B = abs(mdsvalue('-\magnetics::Btor')); % toroidal B-field [T], abs() to keep positive for reversed field
    t_B = mdsvalue('dim_of(-\magnetics::Btor)'); % time
    
    if strcmp(B(1), 'J') || strcmp(t_B(1), 'J') || length(B) == 1 || length(t_B) == 1
        % if B,t_B are not callable or only one point
        return 
    end
    
    B0 = interp1(t_B, B, timebase); % [T], interpolated B

mdsclose;

mdsopen('electrons', shot);
    
    NL_04 = mdsvalue('\ELECTRONS::TOP.TCI.RESULTS:NL_04'); % line-averaged density [m^-3]
    t_04 = mdsvalue('dim_of(\ELECTRONS::TOP.TCI.RESULTS:NL_04)'); % time for density
    
    if strcmp(NL_04(1), 'J') || strcmp(t_04(1), 'J') || length(NL_04) == 1 || length(t_04) == 1
        % if NL_04,t_04 are not callable or only one point
        return 
    end
    
    n20 = 1e-20*interp1(t_04, NL_04, timebase); % interpolated line-averaged density [1e20 m^-3]

mdsclose;

I_M = get_Ip(shot, timebase)*1e-6; %[MA]
k = get_elongation(shot, timebase); % elongation

[~,~,~,~,~, p_input] = get_power(shot, timebase); % p_input = p_ohm + p_lh + p_rf
P_M = p_input*1e-6; % [MW] input power in MW

tE = Wp./(p_input); % energy confinement time (s)

A = get_avg_isotopic_mass_number(shot, timebase); % average isotopic mass number, default = 2

% ITER-89 L-Mode scaling law for confinement time from the Magnetic Fusion
% Energy Formulary (Hartwig)
t89 = 0.048.*(I_M.^0.85).*(R0.^1.2).*(a.^0.3).*(k.^0.5).*(n20.^0.1).*(B0.^0.2).*(A.^0.5).*(P_M.^-0.5);

H89 = real(tE./t89); % get real part, often times this is imaginary, but the imaginary components are all 0 (?)

% reorient all output arrays to original configuration
H89 = orient_array(H89, original_orientation);

end