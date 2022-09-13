function Te_HWHM = get_TS_data_cmod(shot, timebase)

% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Output:
%   Te_HWHM(timebase) = half-width at half-max of parabolic fits to the
%     Te(t,z) data (in meters)
%
% Authors: Robert Granetz, Cristina Rea, & Kevin Montes        May 2017

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
else;
  timebase_column = transpose(timebase);
end;

[~, status] = mdsopen('electrons',shot);
if (mod(status,2)==0)
    Te_HWHM = NaN(size(timebase));
    return
end

% Read in Thomson core temperature data, which is a 2-D array, with the
% dependent dimensions being time and z (vertical coordinate)

node = '.yag_new.results.profiles';

[~,statA] = mdsvalue([node ':te_rz']);
[~,statB] = mdsvalue(['dim_of(' node ':te_rz)']);
[~,statC] = mdsvalue([node ':z_sorted']);

if (mod(statA,2) == 1 & mod(statB,2) == 1 & mod(statC,2) == 1)
    TS_data = mdsvalue([node ':te_rz']);
    TS_time = mdsvalue(['dim_of(' node ':te_rz)']);
    TS_z = mdsvalue([node ':z_sorted']);
    mdsclose;
    
    Te_HWHM = NaN(length(TS_time),1);
    
    % Now, compute Te profile width at each time slice.
    
    zarray = [-0.5:.001:0.5];
    itimes = find(TS_time > 0);
    
    for i = 1:length(itimes);
        y = TS_data(itimes(i),:);
        ok_indices = find(y~=0);
        y = y(ok_indices);
        z = TS_z(ok_indices);
        if (length(ok_indices) > 2)
            %
            [sigma,~,~] = gaussian_fit(z,y);
            Te_HWHM(itimes(i)) = sigma*1.1774;
            %{
            % Code for parabolic fit
            p = polyfit(z, transpose(y), 2);
            if p(1) < 0 % Sanity check -> make sure fit is concave down
                Te_array = polyval(p, zarray);
                [Te_max, maxindx] = max(Te_array);
                z_max = zarray(maxindx); % Location of T profile peak
                Te_HM = Te_max/2;
                [~, HM_indices] = min(abs(Te_array - Te_HM));
                HM_indx = max(HM_indices);
                z_HM = zarray(HM_indx); % Location at half max of T profile    
                Te_HWHM(itimes(i)) = abs(z_HM - z_max);
            end
            %}
        end;
    end;
    Te_HWHM = interp1(TS_time, Te_HWHM, timebase_column, 'linear');
else;
  mdsclose;
  Te_HWHM = NaN(size(timebase));
  return;
end;

% The output signals are currently all column vectors.  However, we desire
% to have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase,2) > 1);
  Te_HWHM = transpose(Te_HWHM);
end;

end