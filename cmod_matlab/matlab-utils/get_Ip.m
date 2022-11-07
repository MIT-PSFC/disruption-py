function [ip1, ip_dot] = get_Ip(shot, timebase)

% This script calculates the derivative of the plasma current, Ip. 
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   ip_dot = dIp/dt, time derivative of plasma current (A/s)
%   ip1 = (interpolated) plasma current (A)
%
% Author: Alex Tinguely 2015-09-14

[timebase, original_orientation] = make_row_array(timebase);

answer = mdsopen('magnetics', shot);

if strcmp(answer, 'Failed')
    ip1 = zeros(1, length(timebase)) + NaN;
    ip_dot = zeros(1, length(timebase)) + NaN;
    return
end

ip = mdsvalue('\magnetics::ip'); % actual/measured Ip, make it positive
t = mdsvalue('dim_of(\magnetics::ip)'); % time for Ip

if strcmp(ip(1), 'J') || strcmp(t(1), 'J') || length(ip) == 1 || length(t) == 1
    % if ip,t are not callable or only one point
    ip1 = zeros(1, length(timebase)) + NaN; % return NaN
    ip_dot = zeros(1, length(timebase)) + NaN;
    return 
end

ip1 = interp1(t, ip, timebase); % interpolate real Ip over given timebase

dip = diff(ip);
dt = diff(t);

t2 = t(1:end-1) + 0.5*dt; % create an array at the midpoints of t (one less value than t)

ipdot = dip./dt;

ip_dot = interp1(t2, ipdot, timebase);

mdsclose;

ip1 = orient_array(ip1, original_orientation);
ip_dot = orient_array(ip_dot, original_orientation);

end