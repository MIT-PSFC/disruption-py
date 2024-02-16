function M = get_avg_isotopic_mass_number(shot, timebase)

% This function calculates the Average Isotopic Mass Number which is used
% in get_H89(), since this value is used in the ITER-89 L-Mode scaling law
% of the confinement time. Because the plasma is dominated by deuterium
% (D) and hydrogen (H), we only consider these two elements. In practice,
% there are other elements injected (whether purposefully or from
% contamination), but we expect them to make a very small contribution to
% the AVERAGE mass number (hopefully).
%
% Here we start with the H/(H+D) ratio, which is the ratio of the number of
% hydrogen atoms to total hydrogen and deuterium atoms. If we let x be the
% ratio, H be the number of H atoms, D be the number of D atoms, and A be
% the average mass number, our derivation is:
%
%       x = H/(H+D) 
%   --> H/D = x/(1-x)
%       M = (1*H + 2*D)/(H+D) since the mass number of H is 1 and of D is 2
%         = (H/D + 2)/(H/D + 1)
%         = [x/(1-x) + 2]/[x/(1-x) + 1]
%         = (x + 2 - 2x)/(x + 1 - x)
%  ---> M = 2-x
%
%  Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   M = average isotropic mass number
%
% Author: Alex Tinguely 2015-09-18

[timebase, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

mdsopen('spectroscopy', shot);

x = mdsvalue('\SPECTROSCOPY::TOP.HALPH_DALPH.ANALYSIS:H_TO_D_RATIO'); % ratio H/(H+D)
t = mdsvalue('dim_of(\SPECTROSCOPY::TOP.HALPH_DALPH.ANALYSIS:H_TO_D_RATIO)'); % time for ratio

if strcmp(x(1), 'J') || strcmp(t(1), 'J') || length(x) == 1 || length(t) == 1
    % if x or t are not callable because of an error, or x,t only has one point
    
    M = 2; % default to all D2, so M = 2
    
else   
    
    M = interp1(t, 2 - x, timebase);
    
end

mdsclose;

% reorient all output arrays to original configuration
M = orient_array(M, original_orientation);

end
