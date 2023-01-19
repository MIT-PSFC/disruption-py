if (exist('time_until_disrupt','var')==0 || exist('time','var')==0);
  fprintf(1,['Call the routine ''retrieve_all_params'' before calling ' ...
             'this routine.\n']);
  return;
end;

indices_no_disrupt = find(isnan(time_until_disrupt));
indices_disrupt = find(~isnan(time_until_disrupt));

shotlist_no_disrupt = unique(shot(indices_no_disrupt));
shotlist_disrupt = unique(shot(indices_disrupt));

duration_no_disrupt = NaN(size(indices_no_disrupt));
duration_disrupt = NaN(size(indices_disrupt));

for i = 1:length(shotlist_no_disrupt);
  duration_no_disrupt(i) = max(time(shot == shotlist_no_disrupt(i)));
end;
for i = 1:length(shotlist_disrupt);
  duration_disrupt(i) = max(time(shot == shotlist_disrupt(i)));
end;

figure( ...
  'Name', 'Histrograms of C-Mod shot lengths', 'NumberTitle','off', ...
  'Menubar', 'none', 'Toolbar', 'none', ...
  'Units', 'normalized', 'Position', [0.10, 0.40, 0.80, 0.55]);
 
bins = [0 : 0.020 : 2.5]; % 10 ms bin width (bin centers specified)
counts = hist(duration_no_disrupt, bins);
axes('Position', [0.10, 0.53, 0.80, 0.37]);
bar(bins, counts, 1, 'green', 'EdgeColor', 'black');
xlim([0,2.5]);
ylimits = ylim;
set(gca, 'Fontsize', 14);
ylabel('# of cases', 'FontSize', 14);
title('Histogram of Alcator C-Mod discharge durations in 2015 campaign', ...
  'Fontsize', 16);
text(0.05, 0.80, 'Non-disruptive', 'color', 'green', ...
  'units', 'normalized', 'fontsize', 16);

hold on;

bins = [0 : 0.010 : 2.5]; % 10 ms bin width (bin centers specified)
counts = hist(duration_disrupt, bins);
axes('Position', [0.10, 0.10, 0.80, 0.37]);
bar(bins, counts, 1.0, 'red', 'EdgeColor', 'black');
h = findobj(gca, 'Type', 'patch');
set(h, 'edgecolor', 'black');
xlim([0,2.5]);
%ylim(ylimits);
ylim([0,50]);
set(gca, 'Fontsize', 14);
xlabel('Discharge duration [s]', 'FontSize', 14);
ylabel('# of cases', 'FontSize', 14);
text(0.05, 0.80, 'Disruptive', 'color', 'red', ...
  'units', 'normalized', 'fontsize', 16);
hold off;

%set(gcf, 'PaperPositionMode', 'manual');
%set(gcf, 'PaperUnits', 'centimeters');
%set(gcf, 'PaperPosition', [1 1 24 12]);
set(gcf, 'PaperPositionMode', 'auto');
print -dpng duration_histograms.png;
