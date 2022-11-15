if (exist('time','var')==0 || exist('time_until_disrupt','var')==0 || ...
    exist('ip','var')==0 || exist('ip_error','var')==0);
  fprintf(1,['Call the routine ''retrieve_all_params'' before calling ' ...
             'this routine.\n']);
  return;
end;

define_indices;
get_ip_sign;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
set(gcf, 'PaperPositionMode', 'auto');

ii = indices_flattop_disrupt_in_flattop;
delta = 5; % 5 kA bin size
bins = [-510 : delta : 110];
tindx = find(time_until_disrupt < 0.0505 & time_until_disrupt > 0.0295);
indx = intersect(ii, tindx);
slice_data = ip_error(indx)/1.e3 .* ipsign(indx);

counts = hist(slice_data, bins);

total_events = sum(counts);
events_below_60_kA = sum(counts(find(bins < -60)));
fraction_below_60_kA = events_below_60_kA / total_events;
fprintf(1,'%i   %i   %5.3f\n', total_events, events_below_60_kA, ...
  fraction_below_60_kA);

bar(bins, counts/sum(counts)/delta, 'hist');
set(gca,'fontsize',12);

xlim([-500, 100]);
ylim([0,.035]);
xlabel('Ip error [kA]', 'FontSize', 14);
ylabel('Probability [kA^{-1}]', 'FontSize', 14);
title('DIII-D: Distribution of Ip error', 'FontSize', 16);

text(0.20, 0.88, '30-50 ms', 'fontsize', 16, ...
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.20, 0.81, 'before', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.20, 0.74, 'disruption', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');

print('d3d_ip_error_disrupt_030-050_hist.png', '-dpng');

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
set(gcf, 'PaperPositionMode', 'auto');

delta = 5; % 5 kA bin size
bins = [-510 : delta : 110];

jj = indices_flattop_no_disrupt;
slice_data = ip_error(jj)/1.e3 .* ipsign(jj);

counts = hist(slice_data, bins);

total_events = sum(counts);
events_below_60_kA = sum(counts(find(bins < -60)));
fraction_below_60_kA = events_below_60_kA / total_events;
fprintf(1,'%i   %i   %5.3f\n', total_events, events_below_60_kA, ...
  fraction_below_60_kA);

bar(bins, counts/sum(counts)/delta, 'hist');
set(gca,'fontsize',12);

xlim([-500, 100]);
ylim([0,.035]);
xlabel('Ip error [kA]', 'FontSize', 14);
ylabel('Probability [kA^{-1}]', 'FontSize', 14);
title('DIII-D: Distribution of Ip error', 'FontSize', 16);

text(0.20, 0.88, 'All', 'fontsize', 16, ...
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.20, 0.81, 'non-disruptive', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.20, 0.74, 'shots', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');

print('d3d_ip_error_no_disrupt_hist.png', '-dpng');
