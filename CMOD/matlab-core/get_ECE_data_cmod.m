function Te_HWHM = get_ECE_data_cmod(shot, timebase)

% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Output:
%   Te_HWHM(timebase) = half-width at half-max of parabolic fits to the
%     Te(t,r) data from the CMod GPC diagnostic 
%
% Authors: Robert Granetz, Cristina Rea, & Kevin Montes        June 2017
%
% The input array, "timebase", can be either be a row vector (1 x n) or a
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

% If electron tree fails to open, return null result

[~, status] = mdsopen('electrons',shot);
if (mod(status,2)==0)
    Te_HWHM = NaN(size(timebase));
    return
end

% Check that at least one channel has data stored. If not, output NaN for
% all values in the Te_HWHM array. 

node = '.ece.gpc_results';
n = 0; m = 0;
for i=1:9
    [t_Te,status1] = mdsvalue(['dim_of(' node '.te:te' num2str(i) ')']);
    [t_r,status2] = mdsvalue(['dim_of(' node '.rad:r' num2str(i) ')']);
    if (mod(status1,2)==1 & mod(status2,2)==1)
        n = length(t_Te);
        m = length(t_r);
        break
    end
end

if (n < 2 | m < 2)
    Te_HWHM = NaN(size(timebase));
    return
end

% Read in GPC temperature data. The GPC diagnostic has 9 channels, each of
% which measure electron temperature for a radius in the plasma that is
% dependent on the magnetic field. Therefore, there are 9 separate data
% nodes in the tree structure where temperature is stored, and 9 additional
% data nodes where each channel's corresponding mapped radius is stored as
% a function of time. Since the time bases for the radial mapping and the
% temperature data are different, I will interpolate the radii onto the GPC
% timebase.

Te_time = NaN(9,n);
Te_data = NaN(9,n);
rad_time = NaN(9,m);
rad_data = NaN(9,m);
radii = NaN(9,n);
indx = 0;

for i=1:9 % for 9 GPC channels
    [Te_data(i,:),status1] = mdsvalue([node '.te:te' num2str(i)]);
    [Te_time(i,:),status2] = mdsvalue(['dim_of(' node '.te:te'...
        num2str(i) ')']);
    if (mod(status1,2)==1 & mod(status2,2)==1)
        indx = i;
        [rad_data(i,:),status1] = mdsvalue([node '.rad:r' num2str(i)]);
        [rad_time(i,:),status2] = mdsvalue(['dim_of(' node ...
            '.rad:r' num2str(i) ')']);
        if (mod(status1,2)==1 & mod(status2,2)==1)
            radii(i,:) = interp1(rad_time(i,:), rad_data(i,:), ... 
                Te_time(i,:),'linear');
        else
            radii(i,:) = NaN(1,n);
        end
    else
        Te_data(i,:) = NaN(1,n); Te_time(i,:) = NaN(1,n);
    end
end

if all(all(isnan(Te_time)))
    Te_HWHM = NaN(size(timebase));
    return
end


% Now, compute Te profile width at each time slice using Gaussian fit
    
Te_HWHM = NaN(1,n);
rarray = [0.5:.001:1];
itimes = find(Te_time(indx,:) > 0);
    
for i = 1:length(itimes);
    y = Te_data(:,itimes(i));
    ok_indices = find(~isnan(y) & y~=0);
    y = y(ok_indices);
    r = radii(ok_indices,itimes(i));
    if (length(ok_indices) > 2)
        %
        [sigma,~,~] = gaussian_fit(r,y);
        Te_HWHM(itimes(i)) = sigma*1.1774;
        %{
        % Code for parabolic fit
        p = polyfit(r, y, 2);
        if p(1) < 0 % Sanity check -> make sure fit is concave down
            Te_array = polyval(p, rarray);
            [Te_max, maxindx] = max(Te_array);
            r_max = rarray(maxindx); % Location of T profile peak
            Te_HM = Te_max/2;
            [~, HM_indices] = min(abs(Te_array - Te_HM));
            HM_indx = max(HM_indices);
            r_HM = rarray(HM_indx); % Location at half max of T profile    
            Te_HWHM(itimes(i)) = abs(r_HM - r_max);
        end
        %}
    end;
end;
Te_HWHM = interp1(Te_time(1,:), Te_HWHM, timebase_column, 'linear');

% The output signals are currently all column vectors.  However, we desire
% to have the output arrays match the shape of "timebase".  Therefore, if
% "timebase" is a row vector, we need to convert all the outputs to row
% vectors.

if (size(timebase,2) > 1);
  Te_HWHM = transpose(Te_HWHM);
end;

end