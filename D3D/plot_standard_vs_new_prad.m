shot = 164716;
mdsconnect('atlas.gat.com');


mdsopen('bolom', shot);
prad_tot = mdsvalue('\prad_tot');
prad_tot_timebase = mdsvalue('dim_of(\prad_tot)')/1.e3;
mdsclose;

smoothing_window = 0.025; % use 25 ms causal smoothing window
a_structure = getbolo(shot, smoothing_window*1.e3);
b_structure = powers(a_structure);
prad_new = b_structure.pwrmix; % watts
prad_new_timebase = a_structure.rawtime; % seconds

%close([1,2]);
figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none', ...
  'PaperPositionMode', 'auto');

plot(prad_new_timebase, prad_new/1.e6, '-b', 'linewidth', 1.5);
xlim([0, 3]);
ylim([0, 5]);
set(gca, 'linewidth', 1.5, 'fontsize', 12)
xlabel('time [s]', 'FontSize', 14);
ylabel('P_{rad} [MW]', 'FontSize', 14);
title('Shot 164716; new P_{rad} vs standard P_{rad}', 'FontSize', 16);

hold on
plot(prad_tot_timebase, prad_tot/1.e6, '-r', 'linewidth', 1.5);

legend({'new prad (25 ms \bf\it{causal}\rm filter)', 'standard prad\_tot'}, ...
  'location', 'northwest');

print('standard_prad_vs_new_prad.png', '-dpng');

smoothing_window = 0.025; % use 25 ms causal smoothing window
a_structure = getbolo(shot, smoothing_window*1.e3);
b_structure = powers(a_structure);
prad_new = b_structure.pwrmix; % watts
prad_new_timebase = a_structure.rawtime; % seconds

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none', ...
  'PaperPositionMode', 'auto');

plot(prad_new_timebase, prad_new/1.e6, '-b', 'linewidth', 1.5);
xlim([2.373, 2.52]);
ylim([-0.5, 70]);
set(gca, 'linewidth', 1.5, 'fontsize', 12)
xlabel('time [s]', 'FontSize', 14);
ylabel('P_{rad} [MW]', 'FontSize', 14);
title('Shot 164716; new P_{rad} vs standard P_{rad}', 'FontSize', 16);

hold on
plot(prad_tot_timebase, prad_tot/1.e6, '-r', 'linewidth', 1.5);
plot([2.465, 2.465], ylim, '-k');

legend({'new prad (25 ms \bf\it{causal}\rm filter)', 'standard prad\_tot'}, ...
  'location', 'northwest');

print('standard_prad_vs_new_prad_expanded.png', '-dpng');
