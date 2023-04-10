function [n_e, dn_dt, Greenwald_fraction] = get_densities(shot, timebase)

% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Output:
%   n_e = average density at center from TCI chord 4 (m ^-3)
%   dn_dt = time derivate of average density n_e (m^-3 s^-1)
%   Greenwald_fraction = n_e/n_G (where n_G is Greenwald density limit)
%
% Authors: Robert Granetz & Kevin Montes               June 2017
% Modified by:
% Ryan Sweeney, 19/17/01
% -added a check for the length of timebase to handle cases where
%  its length is 1, and interpolation fails.
% -added a check on the length of t_n and n_e (Thomson time and density)
%  to make sure they are the same length.  						 
%
% The input array, "timebase", can be either a row vector (1 x n) or a
% column vector (n x 1).  We desire to have the output arrays match the
% shape of "timebase".  However, Matlab calls to the routine "mdsvalue"
% return column vectors for 1-D signals, so it is simpler to work with
% column vectors within this routine, and then, if necessary, convert the
% outputs to row vectors just before exiting this routine.  So the first
% step is to create a copy of "timebase" that is guaranteed to be a column
% vector.


  
if (size(timebase,1) > 1);
  timebase_column = timebase;
else
  timebase_column = transpose(timebase);
end;

n_e = NaN(size(timebase_column));
dn_dt = NaN(size(timebase_column));
Greenwald_fraction = NaN(size(timebase_column));

% added the following check to exit when the timebase is too
% short to interpolate
if length(timebase) == 1;
    return
end

% Initialize all output arrays to NaN (Not-a-Number) column vectors




% Now read in the average density, n_e, and its time derivative dn_dt
% from the electrons tree. 

[~,status] = mdsopen('electrons',shot);
if (mod(status,2) == 1)
    n_e = mdsvalue('.tci.results:nl_04/0.6');
    [t_n, status] = mdsvalue('dim_of(.tci.results:nl_04)');

    if length(n_e) ~= length(t_n);
        n_e = NaN(size(timebase_column));
        dn_dt = NaN(size(timebase_column));
        Greenwald_fraction = NaN(size(timebase_column));
        return;
    end

    if (mod(status,2) == 1)
        dn_dt = gradient(n_e, t_n);
        n_e = interp1(t_n, n_e, timebase_column, 'linear');
        dn_dt = interp1(t_n, dn_dt, timebase_column, 'linear');
    end  
end
mdsclose; 

% Next, read in the measured plasma current, Ip, from the magnetics tree.
% Then, read in the plasma minor radius, a, from the analysis tree, and
% use both to calculate the Greenwald density at each time step in 
% timebase_column. 

[~, status] = mdsopen('magnetics',shot);
if (mod(status,2) == 1)
  ip = mdsvalue('\ip');
  ip = -ip/1e6; % Convert from A to MA and take positive value
  [t_ip, status] = mdsvalue('dim_of(\ip)');
  if (mod(status,2) ==1)
      mdsclose;
      [~, status] = mdsopen('analysis',shot);
      if (mod(status,2) ==1)
          a = mdsvalue('.efit.results.a_eqdsk:aminor');
          [t_a, status] = mdsvalue('dim_of(.efit.results.a_eqdsk:aminor)');
          if (mod(status,2) ==1)
              ip = interp1(t_ip, ip, timebase_column, 'linear');
              a = interp1(t_a, a, timebase_column, 'linear');
              n_G = ip./(pi*a.^2)*1e20; % Greenwald density in m ^-3
              Greenwald_fraction = abs(n_e./n_G);
          end
      end
      mdsclose;
  else
      mdsclose;
  end
end
