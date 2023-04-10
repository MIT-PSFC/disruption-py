function [ip_error, ip_prog, ip] = get_Ip_error(shot, timebase)

% This script calculates Ip_error = (Ip_programmed - Ip),
% or how much the actual plasma current differs from the requested current.
% It linearly interpolates both the programmed and measured plasma currents
% over the given timebase.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   ip_error = Ip_error = (Ip_programmed - Ip) (A)
%   ip_prog = programmed or requested plasma current (A)
%   ip = actual or measured plasma current (A)
%
% NOTE: Since we will also be storing Ip, we can later normalize this error
% to the plasma current, but for now we will just calculate the difference.
%
% Author: Alex Tinguely 2015-09-25

[timebase, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

mdsopen('pcs',shot);

array = []; % used to record with segments in the pcs are turned "on"

ip_error = zeros(1, length(timebase)) + NaN; % Initialize ip_error to NaN
%ip_prog = zeros(1, length(timebase)) + NaN; % Initialize ip_prog to NaN
%ip = zeros(1, length(timebase)) + NaN; % Initialize ip to NaN
ip = get_Ip(shot, timebase);

for i = 1:4 % for each segment, i.e. seg_01, seg_02, seg_03, seg_04
    
    state = mdsvalue(['getnci(".seg_0', num2str(i), ':seg_num","state")']);
    % The result will be "0" for node "ON" and "1" for node "OFF".
    
    if state == 0
        
       array = [array,i]; % adds segment number to array if node is "ON"
       
    end
    
end

l = length(array);

% Initialize arrays
start_time_array = zeros(1,l); % array for start times
index_array = zeros(1,l); % array for indeces of the times in timebase closest to the start times
ip_interp_array = cell(1,l); % cell array for the interpolated ip_programmed
ip_prog = []; % programmed plasma current for the given timebase and each segment

for i = 1:l % for each segment that is turned "ON"
    
    wire = 0; % initialize
    
    for j = 1:16 % for each of the 16 wires

        name = mdsvalue(['.seg_0', num2str(array(i)), ':p_', num2str(17-j,'%02i'), ':name']);
        % provided name of segment (working backwards from 16, since 16 is usually IP)
        % IMPORTANT NOTE: even if the segment is ON, and there IS a
        % programmed plasma current, there may be ZERO effective plasma
        % current as determined by pid_gains (below). In the PCS GUI, these
        % are distinguised as follows: the actual programmed plasma current is given by
        % "IP", but the zero effective programmed plasma current is given
        % "ip". HOWEVER, the tree does not distinguish these (case
        % insensitive) and both are set to capital "IP". That is why we
        % need to look at pid_gains to determine this.
        
        if strcmp(name, 'IP') % if this wire is plasma current, IP
            
            wire = num2str(17-j,'%02i'); % record the wire
            
            pid_gains = mdsvalue(['.seg_0', num2str(array(i)), ':p_', wire, ':pid_gains']);
            % "proportional", "integral", "derivative" gains of IP signal, 
            % if all zero, indicates that programmed IP is NOT actually
            % activated, so "ip"

            inds = find(pid_gains > 0, 1); % first index where matrix elements of pid_gains is greater than zero
            % if inds is NOT empty, then this is a real programmed IP. 
            % if inds is empty, then  we will set these values to NaN
            
            time = mdsvalue(['dim_of(.seg_0', num2str(array(i)), ':p_', wire,')']); % time for IP
            
            if isempty(inds) % if all pid_gains are zero, so zero "ip"
               
                ip_pcs = zeros(1, length(time)) + NaN; % set this ip_pcs to all NaN
                
            else % if the pid_gains are not zero, so real "IP"
                
                ip_pcs = mdsvalue(['.seg_0', num2str(array(i)), ':p_', wire]); % programmed ip
                
            end
            
            break
            
        end
        
    end
    
    if wire == 0 % if no wire is IP
        
        disp(['Cannot find IP in Segment 0', num2str(i)])
        %ip_error = zeros(1, length(timebase)) + NaN; % make ip_error all NaN
        continue % onto the next ACTIVE segment
    
    end    
    
    start_time_array(i) = mdsvalue(['.seg_0', num2str(array(i)), ':start_time']); % start time for this segment

    [~, index_array(i)] = min(abs(start_time_array(i) - timebase)); % index of time in timebase closest to the start time

    ip_interp_array{i} = interp1(time, ip_pcs, timebase); % ip_programmed interpolated over given timebase
                 
end

if  ~exist('wire', 'var') % if the variable 'wire' doesn't even exist (likely from a bad shot)
    return
end

if wire == 0 % after going through all of the active segments, if no wire for IP
             % is found, then we return the array of NaNs
    return
end

[~, indeces] = sort(start_time_array); % put the start times in order (least to greatest)
% and get the indeces of the segments in that order
    
for i = 1:l-1 % for the earliest segment to the second to last segment
    
    % add the subset of the interpolated programmed ip between its start
    % time and the following segment's start time (minus 1 index)
    ip_prog = [ip_prog, ip_interp_array{indeces(i)}(index_array(indeces(i)):index_array(indeces(i+1))-1)];
    
end

% for the last segment, go from its start time to the end of the array
ip_prog = [ip_prog, ip_interp_array{indeces(l)}(index_array(indeces(l)):end)];

mdsclose;

ip_error = ip_prog - ip; % calculate the ip_error

% reorient all output arrays to original configuration
ip = orient_array(ip, original_orientation);
ip_prog = orient_array(ip_prog, original_orientation);
ip_error = orient_array(ip_error, original_orientation);

end