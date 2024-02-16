if (exist('radiated_fraction','var')==0 || ...
    exist('time_until_disrupt','var')==0);
  fprintf(1,['Call the routine ''retrieve_all_params'' before calling ' ...
             'this routine.\n']);
  return;
end;

get_flattop_indices;

no_disrupt_flattop = ...
  intersect(indices_flattop, find(isnan(time_until_disrupt)));

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
plot(time(no_disrupt_flattop), radiated_fraction(no_disrupt_flattop),'.');
xlim([0,2.0]);
ylim([0,2]);
set(gca, 'fontsize', 12);
xlabel('Time [s]', 'FontSize', 14);
ylabel('Radiated power fraction', 'FontSize', 14);
title('Radiated power fraction for non-disruptive shots (flattop ONLY)', ...
  'FontSize', 16);
 
set(gcf, 'PaperPositionMode', 'auto');
print -dpng pradfrac_plot_1.png;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
plot(-time_until_disrupt(indices_flattop), ...         
  radiated_fraction(indices_flattop),'.');
xlim([-1,0]);
ylim([0,2]);
set(gca, 'fontsize', 12);
xlabel('Time before disrupt [s]', 'FontSize', 14);
ylabel('Radiated power fraction', 'FontSize', 14);
title('Radiated power fraction vs time before disrupt (flattop ONLY)', 'FontSize', 16);
 
set(gcf, 'PaperPositionMode', 'auto');
print -dpng pradfrac_plot_2.png;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
plot(-smear(time_until_disrupt(indices_flattop))*1.e3, ...         
  radiated_fraction(indices_flattop),'.');
xlim([-50,0]);
ylim([0,2]);
set(gca, 'fontsize', 12);
xlabel('Time before disrupt [ms]', 'FontSize', 14);
ylabel('Radiated power fraction', 'FontSize', 14);
title('Radiated power fraction vs time before disrupt (flattop ONLY)', 'FontSize', 16);
 
set(gcf, 'PaperPositionMode', 'auto');
print -dpng pradfrac_plot_3.png;

bins = -0.1 : 0.05 : 2.1 ;
counts = hist(radiated_fraction(no_disrupt_flattop), bins);
figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
bar(bins, counts/sum(counts));
xlim([0,2]);
ylim([0,0.25]);
set(gca, 'fontsize', 12);
xlabel('Radiated power fraction', 'FontSize', 14);
ylabel('Probability', 'FontSize', 14);
title(['Histogram of radiated power fraction, ' ...
  'non-disruptive, flattop ONLY'], 'FontSize', 16);

set(gcf, 'PaperPositionMode', 'auto');
print -dpng pradfrac_plot_4.png;

indices = intersect(indices_flattop, find(time_until_disrupt > 0.020));
counts = hist(radiated_fraction(indices), bins);
figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
bar(bins, counts/sum(counts));
xlim([0,2]);
ylim([0,0.25]);
set(gca, 'fontsize', 12);
xlabel('Radiated power fraction', 'FontSize', 14);
ylabel('Probability', 'FontSize', 14);
title(['Histogram of radiated power fraction, ' ...
  'more than 20 ms before disruption, flattop ONLY'], 'FontSize', 16);

set(gcf, 'PaperPositionMode', 'auto');
print -dpng pradfrac_plot_5.png;

indices = intersect(indices_flattop, find(time_until_disrupt <= 0.020));
counts = hist(radiated_fraction(indices), bins);
figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
bar(bins, counts/sum(counts));
xlim([0,2]);
ylim([0,0.25]);
set(gca, 'fontsize', 12);
xlabel('Radiated power fraction', 'FontSize', 14);
ylabel('Probability', 'FontSize', 14);
title(['Histogram of radiated power fraction, ' ...
  'within 20 ms of disruption, flattop ONLY'], 'FontSize', 16);

set(gcf, 'PaperPositionMode', 'auto');
print -dpng pradfrac_plot_6.png;

