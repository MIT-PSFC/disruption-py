if (exist('time','var')==0 || exist('time_until_disrupt','var')==0 || ...
    exist('ip','var')==0 || exist('ip_error','var')==0);
  fprintf(1,['Call the routine ''retrieve_all_params'' before calling ' ...
             'this routine.\n']);
  return;
end;

get_ip_sign;
get_flattop_indices;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
bins = [-700 : 5 : 200]; % 5 kA bin size
tindx = find(time_until_disrupt < 0.0155 & time_until_disrupt > 0.0095);
indx = intersect(indices_flattop, tindx);
slice_data = ip_error(indx)/1.e3 .* ipsign(indx);

counts = hist(slice_data, bins);

total_events = sum(counts);
events_below_70_kA = sum(counts(find(bins < -70)));
fraction_below_70_kA = events_below_70_kA / total_events;
fprintf(1,'%i   %i   %5.3f\n', total_events, events_below_70_kA, ...
  fraction_below_70_kA);

bar(bins, counts/sum(counts)/5, 'hist');
set(gca,'fontsize',11);

xlim([-350, 100]);
ylim([0, 0.03]);
xlabel('Ip error [kA]', 'FontSize', 14);
ylabel('Probability [kA^{-1}]', 'FontSize', 14);
title('Distribution of Ip error (flattop ONLY)', 'FontSize', 16);

text(0.20, 0.87, '[10-15] ms', 'fontsize', 16, ...
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.20, 0.81, 'before', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.20, 0.75, 'disruption', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');

set(gcf, 'PaperPositionMode', 'auto');
%print -dpng ip_error_hist_1.png;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
bins = [-700 : 5 : 200]; % 5 kA bin size

indices_no_disrupt = find(~isnan(time_until_disrupt));
indx = intersect(indices_flattop, indices_no_disrupt);
slice_data = ip_error(indx)/1.e3 .* ipsign(indx);

counts = hist(slice_data, bins);

total_events = sum(counts);
events_below_70_kA = sum(counts(find(bins < -70)));
fraction_below_70_kA = events_below_70_kA / total_events;
fprintf(1,'%i   %i   %5.3f\n', total_events, events_below_70_kA, ...
  fraction_below_70_kA);

bar(bins, counts/sum(counts)/5, 'hist');
set(gca,'fontsize',11);

xlim([-350, 100]);
ylim([0, 0.03]);
xlabel('Ip error [kA]', 'FontSize', 14);
ylabel('Probability [kA^{-1}]', 'FontSize', 14);
title('Distribution of Ip error (flattop ONLY)', 'FontSize', 16);

text(0.20, 0.87, 'All', 'fontsize', 16, ...
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.20, 0.81, 'non-disruptive', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.20, 0.75, 'shots', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');

set(gcf, 'PaperPositionMode', 'auto');
%print -dpng ip_error_hist_2.png;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
bins = [-700 : 5 : 200]; % 5 kA bin size

indices_no_disrupt = find(~isnan(time_until_disrupt));
indx = intersect(indices_flattop, indices_no_disrupt);
indx = intersect(indx, find(time > 0.6));
slice_data = ip_error(indx)/1.e3 .* ipsign(indx);

counts = hist(slice_data, bins);

total_events = sum(counts);
events_below_70_kA = sum(counts(find(bins < -70)));
fraction_below_70_kA = events_below_70_kA / total_events;
fprintf(1,'%i   %i   %5.3f\n', total_events, events_below_70_kA, ...
  fraction_below_70_kA);

bar(bins, counts/sum(counts)/5, 'hist');
set(gca,'fontsize',11);

xlim([-350, 100]);
ylim([0, 0.03]);
xlabel('Ip error [kA]', 'FontSize', 14);
ylabel('Probability [kA^{-1}]', 'FontSize', 14);
title('Distribution of Ip error (flattop ONLY)', 'FontSize', 16);

text(0.20, 0.87, 'All', 'fontsize', 16, ...
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.20, 0.81, 'non-disruptive', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.20, 0.75, 'shots', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');

set(gcf, 'PaperPositionMode', 'auto');
%print -dpng ip_error_hist_2.png;
