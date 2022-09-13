function plot_data(data)

% This function plots data from the 'data' cell array output by the routine
% velocity_data.m. 
%
% Author: Kevin Montes
% Date: 04/14/2017
%
%-------------------------------------------------------------------------%
t = cell2mat(data(:,2));
t_until_disrupt = cell2mat(data(:,3));
var = cell2mat(data(:,4));
var_name = data{1,5};


figure; xlabel('Time(s)'); ylabel([var_name, ' (Default Units)']);
hold
for i=1:size(data,1)
    time = cell2mat(data(i,2));
    v = cell2mat(data(i,4));
    if all(isnan(cell2mat(data(i,3))))
        plot(time,v,'r')
    else
        plot(time,v,'b')
    end
end

h = zeros(2,1);

h(1) = plot(NaN,NaN,'r');
h(2) = plot(NaN,NaN,'b');
legend(h,'No Disruption','Disrupted','Location','best')

figure; 
plot(-t_until_disrupt,var); 
xlabel('Time Until Disruption (s)')
ylabel([var_name, ' (Default Units)'])