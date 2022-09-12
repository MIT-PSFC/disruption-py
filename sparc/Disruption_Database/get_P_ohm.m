function [P_ohm] = get_P_ohm(timebase, original_orientation, smooth_span)

% NOTE: While this function still works, Bob Granetz has written a better
% version in his directory. That one will be used for the disruption
% warning database, so I will leave this unchanged currently.
% Alex Tinguely, 2015-10-13

% This function calculates the ohmic power, P_ohm, in the plasma. We have
% two cases: For Case 1, the average spacing between times in timebase is
% greater than 5ms (EFIT times are every 20ms). Then we can just use
% normal EFIT data for this calculation.
% This is adapted from a scope (/home/tinguely/scopes/shot_basics.dat):
%
%   _vsurf =  deriv(smooth1d(\ANALYSIS::EFIT_SSIBRY,2))*$2pi;
%   _ip=abs(\ANALYSIS::EFIT_AEQDSK:CPASMA);
%   _li = \analysis::efit_aeqdsk:ali;
%   _L = _li*6.28*67.*1.e-9;
%   _vi = _L*deriv(smooth1d(_ip,2));
%   _poh=_ip*(_vsurf-_vi)/1.e6 [MW]
%
% If no smooth_span argument is included, I set the default to 1, which is no
% smoothing. As seen above, this is normally set to 2, and I'm not sure how
% the smooth1d() function works. For Matlab, the span needs to be ODD, so 2
% would automatically be set to 1. I prefer no smoothing since we are
% already on EFIT times at 50 Hz (20ms intervals) which is pretty long. 
%
% For Case 2, the average spacing between times in timebase is less than
% 5ms, so we need a finer time resolution for our data. We use EFIT01,
% which measures at 1kHz, so 1ms time intervals. This calculation comes
% from Ian Faust (very similar to the one above):
%
%   _x=50;_vsurf = deriv(smooth1d(\EFIT01::EFIT_SSIBRY,_x))*$2pi;
%   _ip=abs(\EFIT01::EFIT_AEQDSK:CPASMA);
%   _li = \EFIT01::efit_aeqdsk:ali;
%   _L = _li*6.28*67.*1.e-9;
%   _vi = _L*deriv(smooth1d(_ip,_x));_poh=_ip*(_vsurf-_vi)/1.e6
%
% If no smooth_span is given for this case, I set the default to 50 ms, as
% recommended by Ian. This should smooth out the 60Hz oscillations. NOTE:
% Sometimes this signal can be VERY noisy and even have a significant time
% jump (20ms or so) in the data, which can lead to some strange results. I
% might always stick with Case 1.
%
% Inputs:
%   timebase = times to calculate P_ohm at 
%   original_orientation = 0 if originally column array, 1 if originally row array
%   smooth_span = OPTIONAL argument, the number of elements to
%        average/smooth over. Default is 1, which means no smoothing. If
%        even, sets to smooth_span-1.
%   
% Outputs:
%   P_ohm = ohmic power [W]
%
% Author: Alex Tinguely 15-09-16

%[timebase, original_orientation] = make_row_array(timebase);

P_ohm = zeros(1, length(timebase)); %initialize to all zeros

%mdsopen('mhd', shot);

if mean(diff(timebase)) > 5e-3 % if the average time interval is greater than 5ms (EFIT is 20ms)
                               % use EFIT data

    if nargin == 3 % if smooth_span is provided
       span = smooth_span;   
    else % if smooth_span is not provided
       span = 1; % no smoothing at all
    end

    psi = mdsvalue('\ANALYSIS::EFIT_SSIBRY'); % flux at the surface, I think
    t_psi = mdsvalue('dim_of(\ANALYSIS::EFIT_SSIBRY)'); % measured at 50Hz for EFIT (20ms intervals)

    % if psi, t_psi are not callable or only one point, return array of
    % zeros
    if strcmp(psi(1), 'J') || strcmp(t_psi(1), 'J') || length(psi) == 1 || length(t_psi) == 1
        return
    end
    
    psi_smooth = smooth(psi, span);

    dpsi = diff(psi_smooth);
    dt_psi = diff(t_psi);
    v_surf = dpsi./dt_psi*2*pi; % surface voltage (not sure where 2 pi comes from)

    t2 = t_psi(1:end-1) + 0.5*dt_psi; % create an array at the midpoints of t_psi (one less value than t_psi)

    i_p = abs(mdsvalue('\ANALYSIS::EFIT_AEQDSK:CPASMA')); % plasma current
    t_i = mdsvalue('dim_of(\ANALYSIS::EFIT_AEQDSK:CPASMA)'); % should be same as t_psi
 
    % if i_p, t_i are not callable or only one point, return zeros
    if strcmp(i_p(1), 'J') || strcmp(t_i(1), 'J') || length(i_p) == 1 || length(t_i) == 1
        return
    end    
    
    ip_smooth = smooth(i_p, span);

    di_p = diff(ip_smooth);
    dt_i = diff(t_i);

    ip_interp = interp1(t_i, i_p, t2); % one less value than i_p

    li = mdsvalue('\analysis::efit_aeqdsk:ali'); % inductance
    t_li = mdsvalue('dim_of(\analysis::efit_aeqdsk:ali)'); % should be same as t_psi

    % if li, t_li are not callable or only one point, return zeros
    if strcmp(li(1), 'J') || strcmp(t_li(1), 'J') || length(li) == 1 || length(t_li) == 1
        return
    end
    
    L = li*6.28*67*1e-9; % = li*2pi*(R = 67 cm)*1e-9
    L_interp = interp1(t_li, L, t2); % one less value than L

    v_i = L_interp.*di_p./dt_i; % voltage = L*dI/dt

    p_ohm = ip_interp.*(v_surf - v_i); % [W], power = current*voltage_difference

    P_ohm = interp1(t2, p_ohm, timebase); % [W]

else % if the time intervals are less than 5ms, use EFIT01 data (1kHz, 1ms)
    
    if nargin == 3 % if smooth_span is provided
       span = smooth_span;   
    else % if smooth_span is not provided
       span = 50; % default to 50ms for EFIT01
    end

    psi = mdsvalue('\EFIT01::EFIT_SSIBRY'); % flux at the surface, I think
    t_psi = mdsvalue('dim_of(\EFIT01::EFIT_SSIBRY)'); % measured at 1Hz for EFIT01 (1ms intervals)

    % if psi, t_psi are not callable or only one point, return zeros
    if strcmp(psi(1), 'J') || strcmp(t_psi(1), 'J') || length(psi) == 1 || length(t_psi) == 1
        return
    end
    
    psi_smooth = smooth(psi, span);

    dpsi = diff(psi_smooth);
    dt_psi = diff(t_psi);
    v_surf = dpsi./dt_psi*2*pi; % surface voltage (not sure where 2 pi comes from)

    t2 = t_psi(1:end-1) + 0.5*dt_psi; % create an array at the midpoints of t_psi (one less value than t_psi)

    i_p = abs(mdsvalue('\EFIT01::EFIT_AEQDSK:CPASMA')); % plasma current
    t_i = mdsvalue('dim_of(\EFIT01::EFIT_AEQDSK:CPASMA)'); % should be same as t_psi

    % if i_p, t_i are not callable or only one point, return zeros
    if strcmp(i_p(1), 'J') || strcmp(t_i(1), 'J') || length(i_p) == 1 || length(t_i) == 1
        return
    end  
    
    ip_smooth = smooth(i_p, span);

    di_p = diff(ip_smooth);
    dt_i = diff(t_i);

    ip_interp = interp1(t_i, i_p, t2); % one less value than i_p

    li = mdsvalue('\EFIT01::efit_aeqdsk:ali'); % inductance
    t_li = mdsvalue('dim_of(\EFIT01::efit_aeqdsk:ali)'); % should be same as t_psi

    % if li, t_li are not callable or only one point, return zeros
    if strcmp(li(1), 'J') || strcmp(t_li(1), 'J') || length(li) == 1 || length(t_li) == 1
        return
    end
    
    L = li*6.28*67*1e-9; % = li*2pi*(R = 67 cm)*1e-9
    L_interp = interp1(t_li, L, t2); % one less value than L

    v_i = L_interp.*di_p./dt_i; % voltage = L*dI/dt

    p_ohm = ip_interp.*(v_surf - v_i); % [W], power = current*voltage_difference

    P_ohm = interp1(t2, p_ohm, timebase); % [W]  
    
end

% reorient all output arrays to original orientation
P_ohm = orient_array(P_ohm, original_orientation);

%mdsclose;

end