if (exist('p_rad','var')==0 || exist('dprad_dt','var')==0);
  retrieve_all_params;
end;

get_flattop_indices;
get_ip_sign;

indx_flattop_disrupt = ...
  intersect(find(~isnan(time_until_disrupt)), indices_flattop);

indx_flattop_non_disrupt = ...
  intersect(find(isnan(time_until_disrupt)), indices_flattop);

indx_shot = find(shot == 1150923014); % from Matt Reinke's e-mail

close all;

figure(1);
plot(-smear(time_until_disrupt(indices_flattop), 0.0005)*1.e3, ...
  p_rad(indices_flattop)/1.e6, '.');
xlim([-40,0]);
ylim([-0.1, 50]);
set(gcf, 'menubar', 'none', 'toolbar', 'none');
set(gca, 'fontsize', 12);
xlabel('Time until disrupt [ms]', 'fontsize', 14);
ylabel('P_{rad} [MW]', 'fontsize', 14);
title('P_{rad} vs time until disrupt', 'fontsize', 16);
text(0.05, 0.90, 'during flattop only', 'color', 'black', ...
  'units', 'normalized', 'fontsize', 16);
set(gcf, 'PaperPositionMode', 'auto');
%print -dpng prad_vs_time_until_disrupt.png;

hold on;
plot(-time_until_disrupt(indx_shot)*1.e3, p_rad(indx_shot)/1.e6, '.r');
hold off;
set(gcf, 'PaperPositionMode', 'auto');
%print -dpng prad_vs_time_until_disrupt_and_1150923014.png;

figure(2);
plot(-smear(time_until_disrupt(indices_flattop), 0.0005)*1.e3, ...
  dprad_dt(indices_flattop)/1.e9, '.');
xlim([-40,0]);
ylim('auto');
set(gcf, 'menubar', 'none', 'toolbar', 'none');
set(gca, 'fontsize', 12);
xlabel('Time until disrupt [ms]', 'fontsize', 14);
ylabel('dP_{rad}/dt [MW/ms]', 'fontsize', 14);
title('dP_{rad}/dt vs time until disrupt', 'fontsize', 16);
text(0.05, 0.90, 'during flattop only', 'color', 'black', ...
  'units', 'normalized', 'fontsize', 16);
set(gcf, 'PaperPositionMode', 'auto');
%print -dpng dprad_dt_vs_time_until_disrupt.png;

figure(3);
plot(-smear(time_until_disrupt(indices_flattop), 0.0005)*1.e3, ...
  dprad_dt(indices_flattop)/1.e9, '.');
xlim([-40,0]);
ylim([-10,10]);
set(gcf, 'menubar', 'none', 'toolbar', 'none');
set(gca, 'fontsize', 12);
xlabel('Time until disrupt [ms]', 'fontsize', 14);
ylabel('dP_{rad}/dt [MW/ms]', 'fontsize', 14);
title('dP_{rad}/dt vs time until disrupt', 'fontsize', 16);
text(0.05, 0.90, 'during flattop only', 'color', 'black', ...
  'units', 'normalized', 'fontsize', 16);
set(gcf, 'PaperPositionMode', 'auto');
%print -dpng dprad_dt_zoomed_vs_time_until_disrupt.png;

hold on;
plot(-time_until_disrupt(indx_shot)*1.e3, dprad_dt(indx_shot)/1.e9, '.r');
hold off;
set(gcf, 'PaperPositionMode', 'auto');
%print -dpng dprad_dt_vs_time_until_disrupt_and_1150923014.png;

figure(4);
binsize = 0.1; % MW
bins = [0 : binsize : 50];
counts = hist(p_rad(indx_flattop_non_disrupt)/1.e6, bins);
bar(bins, counts/sum(counts)/binsize, 'hist');
xlim([0, 50]);
ylim('auto');
set(gcf, 'menubar', 'none', 'toolbar', 'none');
set(gca, 'fontsize', 12);
xlabel('P_{rad} [MW]', 'fontsize', 14);
ylabel('Probability [MW]^{-1}', 'fontsize', 14);
title({'Distribution of P_{rad}', 'for non-disruptions'}, 'fontsize', 16);
text(0.05, 0.90, 'during flattop only', 'color', 'black', ...
  'units', 'normalized', 'fontsize', 16);
set(gcf, 'PaperPositionMode', 'auto');
%print -dpng histogram_prad_non_disrupt.png;

figure(5);
binsize = 0.1; % MW
bins = [0 : binsize : 50];
counts = hist(p_rad(indx_flattop_non_disrupt)/1.e6, bins);
bar(bins, counts/sum(counts)/binsize, 'hist');
xlim([0, 50]);
ylim([0, 0.01]);
set(gcf, 'menubar', 'none', 'toolbar', 'none');
set(gca, 'fontsize', 12);
xlabel('P_{rad} [MW]', 'fontsize', 14);
ylabel('Probability [MW]^{-1}', 'fontsize', 14);
title({'Distribution of P_{rad}', 'for non-disruptions'}, 'fontsize', 16);
text(0.05, 0.90, 'during flattop only', 'color', 'black', ...
  'units', 'normalized', 'fontsize', 16);
set(gcf, 'PaperPositionMode', 'auto');
%print -dpng histogram_prad_non_disrupt_zoomed.png;

figure(6);
binsize = 0.1; % MW
bins = [0 : binsize : 50];
counts = hist(p_rad(indx_flattop_disrupt)/1.e6, bins);
bar(bins, counts/sum(counts)/binsize, 'hist');
xlim([0, 50]);
ylim([0, 0.01]);
set(gcf, 'menubar', 'none', 'toolbar', 'none');
set(gca, 'fontsize', 12);
xlabel('P_{rad} [MW]', 'fontsize', 14);
ylabel('Probability [MW]^{-1}', 'fontsize', 14);
title({'Distribution of P_{rad}', 'for flattop disruptions'}, ...
  'fontsize', 16);
text(0.05, 0.90, 'during flattop only', 'color', 'black', ...
  'units', 'normalized', 'fontsize', 16);
set(gcf, 'PaperPositionMode', 'auto');
%print -dpng histogram_prad_disrupt.png;

figure(7);
binsize = 0.1; % MW/ms
bins = [-80 : binsize : 80];
counts = hist(dprad_dt(indx_flattop_non_disrupt)/1.e9, bins);
bar(bins, counts/sum(counts)/binsize, 'hist');
xlim([-10,10]);
ylim('auto');
set(gcf, 'menubar', 'none', 'toolbar', 'none');
set(gca, 'fontsize', 12);
xlabel('dP_{rad}/dt [MW/ms]', 'fontsize', 14);
ylabel('Probability [MW/ms]^{-1}', 'fontsize', 14);
title({'Distribution of dP_{rad}/dt', 'for non-disruptions'}, 'fontsize', 16);
text(0.05, 0.90, 'during flattop only', 'color', 'black', ...
  'units', 'normalized', 'fontsize', 16);
set(gcf, 'PaperPositionMode', 'auto');
%print -dpng histogram_dprad_dt_non_disrupt.png;

figure(8);
binsize = 0.1; % MW/ms
bins = [-80 : binsize : 80];
counts = hist(dprad_dt(indx_flattop_non_disrupt)/1.e9, bins);
bar(bins, counts/sum(counts)/binsize, 'hist');
xlim([-10,10]);
ylim([0, 0.01]);
set(gcf, 'menubar', 'none', 'toolbar', 'none');
set(gca, 'fontsize', 12);
xlabel('dP_{rad}/dt [MW/ms]', 'fontsize', 14);
ylabel('Probability [MW/ms]^{-1}', 'fontsize', 14);
title({'Distribution of dP_{rad}/dt', 'for non-disruptions'}, 'fontsize', 16);
text(0.05, 0.90, 'during flattop only', 'color', 'black', ...
  'units', 'normalized', 'fontsize', 16);
set(gcf, 'PaperPositionMode', 'auto');
%print -dpng histogram_dprad_dt_non_disrupt_zoomed.png;

figure(9);
binsize = 0.1; % MW/ms
bins = [-80 : binsize : 80];
counts = hist(dprad_dt(indx_flattop_disrupt)/1.e9, bins);
bar(bins, counts/sum(counts)/binsize, 'hist');
xlim([-10,10]);
ylim([0, 0.01]);
set(gcf, 'menubar', 'none', 'toolbar', 'none');
set(gca, 'fontsize', 12);
xlabel('dP_{rad}/dt [MW/ms]', 'fontsize', 14);
ylabel('Probability [MW/ms]^{-1}', 'fontsize', 14);
title({'Distribution of dP_{rad}/dt', 'for flattop disruptions'}, ...
  'fontsize', 16);
text(0.05, 0.90, 'during flattop only', 'color', 'black', ...
  'units', 'normalized', 'fontsize', 16);
set(gcf, 'PaperPositionMode', 'auto');
%print -dpng histogram_dprad_dt_disrupt.png;
