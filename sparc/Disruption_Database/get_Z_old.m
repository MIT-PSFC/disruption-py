function [Z_cur1, vZ, Z_times_vZ] = get_Z_old(shot, timebase)

% This script calculates the vertical velocity, vZ, from the vertical
% position, Z_cur. Finally, this function also outputs the product of Z and vZ,
% which is large and positive when the plasma is moving upwards above the
% midplane and moving downwards below the midplane.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   vZ = vertical velocity (m/s)
%   Z_times_vZ = Z_cur1*vZ (m^2/s)
%   Z_cur1 = interpolated vertical position (m)
%
% Author: Alex Tinguely 2015-09-09

[timebase_row, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

[~,status] = mdsopen('mhd', shot);

if mod(status, 2) == 0 % if this shot does not exist, there is an error
    
    Z_cur1 = NaN(size(timebase));
    vZ = NaN(size(timebase));
    Z_times_vZ = NaN(size(timebase));
    return    
    
else

    [Z_cur, status1] = mdsvalue('\analysis::efit_aeqdsk:zcur'); % actual/measured ZCUR from EFIT
    [t, status2] = mdsvalue('dim_of(\analysis::efit_aeqdsk:zcur)'); % time for ZCUR

    if mod(status1, 2) == 0 ||mod(status2, 2) == 0  % if an error occurs

        %disp('Error in calling Z_cur');
        Z_cur1 = NaN(size(timebase)); % set Z_cur to NaN
        vZ = NaN(size(timebase)); % set vZ to NaN
        Z_times_vZ = NaN(size(timebase)); % set Z*vZ to NaN
        return

    end 

    Z_cur1 = interp1(t, Z_cur, timebase_row); % interpolate ZCUR over given timebase
    % NOTE: We probably don't have to interpolate because the requested
    % timebase is likely the same as the EFIT timebase, but just in case we
    % choose something different, we have the option.

    dZ = diff(Z_cur);
    dt = diff(t);

    vz = dZ./dt;

    t2 = t(1:end-1) + 0.5*dt; % create an array at the midpoints of t (one less value than t)

    vZ = interp1(t2, vz, timebase_row);

    Z_times_vZ = Z_cur1.*vZ; 

end

mdsclose;

% reorient all output arrays to original configuration
Z_cur1 = orient_array(Z_cur1, original_orientation);
vZ = orient_array(vZ, original_orientation);
Z_times_vZ = orient_array(Z_times_vZ, original_orientation);

end