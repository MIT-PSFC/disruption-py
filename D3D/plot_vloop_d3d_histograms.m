if (exist('time','var')==0 || exist('time_until_disrupt','var')==0 || ...
    exist('v_loop_efit','var')==0);
  fprintf(1,['Call the routine ''retrieve_all_params'' before calling ' ...
             'this routine.\n']);
  return;
end;

define_indices;
get_ip_sign;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
set(gcf, 'PaperPositionMode', 'auto');

delta = 0.2; %0.2 volt bin size
bins = [-6.5 : delta : 6.5];
tindx = find(time_until_disrupt > 0.029 & time_until_disrupt < 0.031);
indx = intersect(indices_flattop_disrupt_in_flattop, tindx);
slice_data = v_loop_efit(indx) .* ipsign(indx);

counts = hist(slice_data, bins);

total_events = sum(counts);
declared_disrupt_events = sum(counts(find(bins < -0.5 | bins > 1.8)));
fraction_declared_disrupt_events = declared_disrupt_events / total_events;
fprintf(1,'%i   %i   %5.3f\n', total_events, declared_disrupt_events, ...
  fraction_declared_disrupt_events);

bar(bins, counts/sum(counts)/delta, 'hist');
set(gca, 'fontsize', 12);

xlim([-6, 6]);
%ylim([0,.035]);
xlabel('V_{loop} from EFIT [volts]', 'FontSize', 14);
ylabel('Probability [volt^{-1}]', 'FontSize', 14);
title('DIII-D: Distribution of V_{loop}', 'FontSize', 16);

text(0.80, 0.91, '30 ms before', 'fontsize', 16, ...
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.80, 0.83, 'flattop', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.80, 0.75, 'disruptions', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');

print('d3d_vloop_disrupt_030_histogram.png', '-dpng');

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
set(gcf, 'PaperPositionMode', 'auto');

delta = 0.2; %0.2 volt bin size
bins = [-6.5 : delta : 6.5];
indx = indices_flattop_no_disrupt;
slice_data = v_loop_efit(indx) .* ipsign(indx);

counts = hist(slice_data, bins);

total_events = sum(counts);
declared_disrupt_events = sum(counts(find(bins < -0.5 | bins > 1.8)));
fraction_declared_disrupt_events = declared_disrupt_events / total_events;
fprintf(1,'%i   %i   %5.3f\n', total_events, declared_disrupt_events, ...
  fraction_declared_disrupt_events);

bar(bins, counts/sum(counts)/delta, 'hist');
set(gca, 'fontsize', 12);

xlim([-6, 6]);
%ylim([0,.035]);
xlabel('V_{loop} from EFIT [volts]', 'FontSize', 14);
ylabel('Probability [volt^{-1}]', 'FontSize', 14);
title('DIII-D: Distribution of V_{loop}', 'FontSize', 16);

text(0.75, 0.91, 'During flattop of', 'fontsize', 16, ...
  'units', 'normalized', 'HorizontalAlignment', 'center');
text(0.75, 0.83, 'non-disruptions', 'fontsize', 16, ... 
  'units', 'normalized', 'HorizontalAlignment', 'center');

print('d3d_vloop_no_disrupt_histogram.png', '-dpng');
