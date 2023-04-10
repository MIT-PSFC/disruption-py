function new_data = num_deriv(data)

%------------------------------------------------------------------------%
% This function takes the cell array 'data' output from velocity_data.m
% and appends an extra column to it populated with the acceleration on
% axis, dv/dt. 
%
% Author: Kevin Montes
% Date Created: 04/14/2017
%
%------------------------------------------------------------------------%

s = size(data,2)+1;
data{1,s} = [];

for i=1:size(data,1)
    t = cell2mat(data(i,2));
    v = cell2mat(data(i,4));
    dvdt = diff(v)./diff(t);
    data{i,s} = [NaN;dvdt];
end

new_data = data;