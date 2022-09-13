function new_data = examine(data,domain,range)

%------------------------------------------------------------------------%
%
% This function filters the 'data' cell array to only include shots with
% velocity values that fall in the specified domain and range. 
%
% Author: Kevin Montes
% Date: 04/14/2017
%
%------------------------------------------------------------------------%
tmin = domain(1);
tmax = domain(2);
vmin = range(1);
vmax = range(2);

indices = 1:size(data,1);

% For time vs. v(t) plot, |tmin|< |tmax|, but for -time_until_disrupt vs. 
% v(t) plot, |tmin|> |tmax|. After assigning t variable to appropriate data
% in loop, can then assign a 0 to each index corresponding to a shot with
% no velocity data falling in the domain and range. 

for i=1:size(data,1)
    if abs(tmin) < abs(tmax)
        t = cell2mat(data(i,2));
    else
        t = -cell2mat(data(i,3));
    end
    v = cell2mat(data(i,4));
    k = find(t<=tmax & t>=tmin);
    if ~any(v(k)>vmin & v(k)<vmax)
        indices(i)=0;
    end
end

new_data = data(indices(indices~=0),:);