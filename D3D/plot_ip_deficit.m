load('ip_pcs_data.mat'); % This loads: ip, ip_pcs, shots, tdis

% The plasma current can be positive or negative, but the pcs programming
% is only positive.  (The buswork has to be switched around to change the
% polarity of the plasma current.)  Therefore the absolute value of Ip has
% to be used when determining the difference from the programmed demand.

ip_deficit = ip_pcs - abs(ip);

set(gcf, 'menubar', 'none');
plot(shots, ip_deficit/1.e3, '.');
xlim([100000, 165000]);
ylim([-800, 1200]);
xtickvals = get(gca, 'XTick');
xticklabels = cell(1,numel(xtickvals));
for itick = 1 : numel(xtickvals);
  xticklabels(itick) = cellstr(sprintf('%i',xtickvals(itick)));
end;
set(gca, 'XTickLabel', xticklabels);
xlabel('Shot number', 'fontsize', 12);
ylabel('Ip deficit (kA)', 'fontsize', 12);
title('Ip(programmed) - Ip(actual) just before disruption', 'fontsize', 14);

figure;
set(gcf, 'menubar', 'none');
xbins = [-1000 : 10 : 1500];
hist(ip_deficit/1.e3, xbins);
xlim([-500, 1000]);
xlabel('Ip(programmed) - Ip(actual) just before disruption [kA]', ...
       'fontsize', 12);
ylabel('Number of occurrences', 'fontsize', 12)';
title('Histogram of Ip deficit just before disruption', 'fontsize', 14);
