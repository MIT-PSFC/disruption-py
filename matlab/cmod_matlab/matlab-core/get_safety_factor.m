function [q0, qstar, q95] = get_safety_factor(timebase, original_orientation)

% This script gets the three safety factors: q0 on axis, qstar = q* for the
% cylindrical approximation, and q95.
%
% Inputs:
%   timebase = array of desired time values
%   original_orientation = 0 if originally column array, 1 if originally row array
%
% Outputs:
%   q0 = safety factor at center
%   qstar = q* = cylindrical safety factor
%   q95 = edge safety factor
%
% Author: Alex Tinguely 2015-09-09

%[timebase, original_orientation] = make_row_array(timebase); % makes timebase into row array
% easier to do calculations this way

%mdsopen('mhd', shot);

[q0, status1] = mdsvalue('\analysis::efit_aeqdsk:q0'); % actual/measured q0
[t0, status2] = mdsvalue('dim_of(\analysis::efit_aeqdsk:q0)'); % time for q0

if mod(status1,2)==0 || mod(status2,2)==0 % if an error occurs
   
    disp('Error in calling q0');
    q0 = NaN(size(timebase)); % set to NaN
     
else
    
    q0 = interp1(t0, q0, timebase); % interpolate measured q0 over given timebase

end

%---------------------------------

[qstar, status1] = mdsvalue('\analysis::efit_aeqdsk:qstar'); % actual/measured q*
[tstar, status2] = mdsvalue('dim_of(\analysis::efit_aeqdsk:qstar)'); % time for q*

if mod(status1,2)==0 || mod(status2,2)==0  % if an error occurs
   
    disp('Error in calling qstar');
    qstar = NaN(size(timebase)); % set to NaN
     
else

    qstar = interp1(tstar, qstar, timebase); % interpolate measured q* over given timebase

end

%----------------------------------

[q95, status1] = mdsvalue('\analysis::efit_aeqdsk:qpsib'); % actual/measured q95
[t95, status2] = mdsvalue('dim_of(\analysis::efit_aeqdsk:qpsib)'); % time for q95

if mod(status1,2)==0 || mod(status2,2)==0 % if an error occurs
   
    disp('Error in calling q95');
    q95 = NaN(size(timebase)); % set to NaN
     
else

    q95 = interp1(t95, q95, timebase); % interpolate q95 over given timebase
    
end

%mdsclose;

% reorient all output arrays to original configuration
q0 = orient_array(q0, original_orientation);
qstar = orient_array(qstar, original_orientation);
q95 = orient_array(q95, original_orientation);

end