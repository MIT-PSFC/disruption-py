function [ne1, ne_dot, nG_ratio] = get_ne(shot, timebase)

% This script calculates the derivative of the plasma density, ne = nebar.
% This also calculates the density normalized to the Greenwald density limit.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   ne1 = (interpolated) density (nebar) (1/m^3))
%   ne_dot = dne/dt, time derivative of the density (1/m^3/s)
%   nG_ratio = ne/nG (m^-3) ratio of density and Greenwald density,
%              Greenwald fraction
%
% Author: Alex Tinguely 2015-09-09

[timebase_row, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

[~, status] = mdsopen('electrons', shot);

if mod(status, 2) == 0 % if this shot does not exist, there is an error
    
    ne1 = NaN(size(timebase));
    ne_dot = NaN(size(timebase));
    nG_ratio = NaN(size(timebase));
    return
    
else

    [ne, status1] = mdsvalue('\electrons::top.tci.results.inversion.nebar_efit'); % actual/measured nebar
    [t, status2] = mdsvalue('dim_of(\electrons::top.tci.results.inversion.nebar_efit)'); % time for nebar

    if mod(status1,2)==0 || mod(status2,2)==0|| length(ne) == 1 || length(t) == 1 
        % if an error occurs or only one point for ne1,t

        disp('Error in calling ne');
        ne1 = NaN(size(timebase));
        ne_dot = NaN(size(timebase));
        nG_ratio = NaN(size(timebase));
        return

    end 

    ne1 = interp1(t, ne, timebase_row); % interpolate real density over given timebase

    dne = diff(ne);
    dt = diff(t);

    nedot = dne./dt;

    t2 = t(1:end-1) + 0.5*dt; % create an array at the midpoints of t (one less value than t)

    ne_dot = interp1(t2, nedot, timebase_row);

    Ip = get_Ip(shot, timebase_row);

    a = 0.22; % minor radius of CMod [m]

    nG = Ip*1e-6/(pi*a^2)*1e20; % Greenwald density limit [m^-3] 
    % Usually, nG is in [1e20 m^-3] and Ip is measured in MA, so I've added the
    % appropriate factors

    nG_ratio = ne1./nG; % density normalized to Greenwald density nG

end

mdsclose;

% reorient all output arrays to original configuration
ne1 = orient_array(ne1, original_orientation);
ne_dot = orient_array(ne_dot, original_orientation);
nG_ratio = orient_array(nG_ratio, original_orientation);

end