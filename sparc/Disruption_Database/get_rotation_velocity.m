function [v_0, v_mid] = get_rotation_velocity(shot, timebase)

% This function returns the rotational velocity of the plasma on the
% magnetic axis (v_0) and at r/a = 0.5 (v_mid). Unfortunately, HIREX does not
% have the velocity at the edge, and while CXRS could theoretically
% calculate it, this has not yet been done, so it is not included in this
% function.
%
% There are several different ways that HIREX_SR determines the rotation
% velocity. We will use these lines:
%
%   A - Lyman alpha from H-like Argon
%   Z - Z line from He-like Argon
%
% According to John Rice, for the velocity on axis, we should use A-line
% ONLY if the counts from the line-integrated intensity are over 500. If
% the count is under 500, then we use the Z-line signal on axis. For the
% velocity at r/a=0.5, we should use Z-line data. If there is absolutely no
% variation in the data (horizontal line, derivative zero), then something
% is wrong. In that case (or in the case of an error), we return an array
% of NaN's. Also, speeds greater than 150 km/s should be removed.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   v_0 = rotation velocity at magnetic axis (r/a = 0) [m/s]
%   v_mid = rotation velocity at r/a = 0.5 [m/s]
%
% Author: Alex Tinguely 2015-09-23

[timebase_row, original_orientation] = make_row_array(timebase);

[~, status] = mdsopen('spectroscopy', shot);

if mod(status,2) == 0 % if we can't open spectroscopy for this shot
    
    v_0 = NaN(size(timebase));
    v_mid = NaN(size(timebase));
    return
    
else

    [v0, status1] = mdsvalue('\SPECTROSCOPY::TOP.HIREX_SR.ANALYSIS:A:VEL'); % [km/s], velocity on axis from A line
                                                                            % this is our default velocity on axis
    v0 = v0*1e3; %[m/s]
    
    [t, status2] = mdsvalue('dim_of(\SPECTROSCOPY::TOP.HIREX_SR.ANALYSIS:A:VEL)'); % time [s]
    [counts_A] = mdsvalue('\SPECTROSCOPY::TOP.HIREX_SR.ANALYSIS:A:INT'); % line-integrated intensity [counts/pixel]

    if mod(status1,2) ==0 || mod(status2,2) ==0 % if there is an error

        counts_A = 0;
        
    end
        
    [vZ, status1Z] = mdsvalue('\SPECTROSCOPY::TOP.HIREX_SR.ANALYSIS:Z:VEL'); % [km/s], velocity MATRIX from Z line
    vZ = vZ*1e3; % [m/s]
    
    vmid = vZ(:,2); % velocity at r/a = 0.5, [m/s]
    [tmid, status2] = mdsvalue('dim_of(\SPECTROSCOPY::TOP.HIREX_SR.ANALYSIS:Z:VEL)'); % time [s], should be same as t

    if mean(counts_A) < 500 % if the average number of counts from A (line-integrated intensity) is less than 500

        v0 = vZ(:,1); % get the velocity on axis from Z-line data, [m/s]
        t = tmid; % change the time to Z time (although likely the same as original t from A-line)
        status1 = status1Z;
        
    end
    
    % initialize the outputs to arrays of NaN if neither is callable
    v_0 = NaN(size(timebase_row));
    v_mid = NaN(size(timebase_row)); 
    
    diff_v0 = max(abs(diff(v0))); % maximum change in v0 (positive)
    diff_vmid = max(abs(diff(vmid))); % maximum change in vmid (positive)

    if ~(mod(status1,2) ==0 || mod(status2,2) ==0  || diff_v0 == 0) 
      % if an error does NOT occur OR if v0 is NOT constant in time

      indeces = find(abs(v0) >= 150e3); 
      % find indeces where the speed is greater than or equal to the 150 km/s = 150,000 m/s (unphysical)

      v0(indeces) = []; % remove those data points
      t(indeces) = [];

      v_0 = interp1(t, v0, timebase_row); % [m/s], interpolate to the timebase (likely the same as the other times)

    end   

    % ------------------- similar calculation for the r/a = 0.5 velocity ------------------------

    if ~(mod(status1Z,2) ==0 || mod(status2,2) ==0 || diff_vmid == 0) 
      % if an error does NOT occur OR if vmid is NOT constant in time

      indeces = find(abs(vmid) >= 150e3); 
      % find indeces where the velocity is greater than or equal to 150,000 m/s

      vmid(indeces) = []; % remove those data points
      tmid(indeces) = [];

      v_mid = interp1(tmid, vmid, timebase_row); % [m/s], interpolate to the timebase (likely the same as the other times)

    end

end

mdsclose;

% reorient all output arrays to original configuration
v_0 = orient_array(v_0, original_orientation);
v_mid = orient_array(v_mid, original_orientation);

end