function [Z_error, Z_prog, Z_cur] = get_Z_error_old(shot, timebase)

% This script calculates Z_error = Z_programmed - Z_cur,
% or how much the actual vertical position differs from the requested position.
% It linearly interpolates both the programmed and measured positions
% over the given timebase. Note that this is just the difference
% so the units are meters (m).
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   Z_error = Z_programmed - Z_cur (m), returns NaN array if no programmed Z
%   Z_prog = programmed or requested Z (m)
%   Z_cur = actual or measured Z (m)
%
% Author: Alex Tinguely 2015-09-09

[timebase_row, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

Z_cur = get_Z(shot, timebase);

array = []; % used to record with segments in the pcs are turned "on"

Z_error = NaN(size(timebase_row)); % Initialize Z_error to NaN

[~, status] = mdsopen('pcs',shot);

if mod(status, 2) == 0 % if this shot does not exist, there is an error
    
    Z_prog = NaN(size(timebase));
    Z_error = NaN(size(timebase));
    return
    
else

    % We need to look in the control system (pcs) at each "segment" to
    % determine where the vertical position (Z) is programmed (if at all).

    for i = 1:4 % for each segment, i.e. seg_01, seg_02, etc.

        state = mdsvalue(['getnci(".seg_0', num2str(i), ':seg_num","state")']);
        % NOTE: The result will be "0" for node "on" and "1" for node "off".

        if state == 0

           array = [array,i]; % adds segment number to array if node is "on"

        end

    end

    l = length(array);

    % Initialize arrays
    start_time_array = zeros(1,l); % array for start times
    index_array = zeros(1,l); % array for indeces of the times in timebase closest to the start times
    Z_interp_array = cell(1,l); % cell array for the interpolated ip_programmed
    Z_prog = []; % programmed plasma current for the given timebase and each segment

    for i = 1:l % for each segment that is turned "on"

        wire = 0; % initialize

        for j = 1:16 % for each of the 16 wires

            name = mdsvalue(['.seg_0', num2str(array(i)), ':p_', num2str(j,'%02i'), ':name']);
            % NOTE: The PCS GUI has different "titles" for ZCUR. When in use
            % there are asterisks around it like "* ZCUR *", and if inactive (I
            % think) it is lower case, "zcur". However, in the Tree, the name
            % is case insensitive, so both active and inactive ZCURs will have
            % the same name "ZCUR". We need to look a the pid_gains to determine
            % whether the effective Z_programmed is zero.

            if strcmp(name, 'ZCUR') % if this wire is vertical position

                wire = num2str(j,'%02i'); % record the wire number

                pid_gains = mdsvalue(['.seg_0', num2str(array(i)), ':p_', wire, ':pid_gains']);
                % "proportional", "integral", "derivative" gains of ZCUR signal, 
                % if all zero, indicates that programmed ZCUR is NOT actually
                % activated, so "zcur"

                inds = find(pid_gains > 0, 1); % first index where matrix elements of pid_gains is greater than zero
                % if inds is NOT empty, then this is a real programmed ZCUR. 
                % if inds is empty, then  we will set these values to NaN

                time = mdsvalue(['dim_of(.seg_0', num2str(array(i)), ':p_', wire,')']); % time

                if isempty(inds) % if all pid_gains are zero

                    Z_pcs = zeros(1, length(time)) + NaN; % set this Z_pcs to all NaN

                else % if the pid_gains are not zero,

                    Z_pcs = abs(mdsvalue(['.seg_0', num2str(array(i)), ':p_', wire])); % programmed ZCUR

                end

                break

            end

        end

        if wire == 0 % if no wire in this active segment is ZCUR

            disp(['Cannot find ZCUR in Segment 0', num2str(i)])
            continue

        end

        start_time_array(i) = mdsvalue(['.seg_0', num2str(array(i)), ':start_time']); % start time for this segment

        [~, index_array(i)] = min(abs(start_time_array(i) - timebase_row)); % index of time in timebase closest to the start time

        Z_interp_array{i} = interp1(time, Z_pcs, timebase_row); % Z_programmed interpolated over given timebase

    end

    if ~exist('wire', 'var') % if no variable 'wire' is defined (usually from a bad shot)
        return               % end the function
    end

    if wire == 0 % after going through all of the active segments, if no wire for ZCUR
                 % is found then we return the array of NaNs
        return
    end

    [~, indeces] = sort(start_time_array); % put the start times in order (least to greatest)
    % and get the indeces of the segments in that order

    for i = 1:l-1 % for the earliest segment to the second to last segment

        % add the subset of the interpolated programmed Z between its start
        % time and the following segment's start time (minus 1 index)
        Z_prog = [Z_prog, Z_interp_array{indeces(i)}(index_array(indeces(i)):index_array(indeces(i+1))-1)];

    end

    % for the last segment, go from its start time to the end of the array
    Z_prog = [Z_prog, Z_interp_array{indeces(l)}(index_array(indeces(l)):end)];

    mdsclose;

end

Z_error = Z_prog - Z_cur; % [m], calculate Z_error

% reorient all output arrays to original configuration
Z_cur = orient_array(Z_cur, original_orientation);
Z_prog = orient_array(Z_prog, original_orientation);
Z_error = orient_array(Z_error, original_orientation);

end