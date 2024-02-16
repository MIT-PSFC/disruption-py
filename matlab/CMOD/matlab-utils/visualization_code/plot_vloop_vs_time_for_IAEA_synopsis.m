db = set_database('logbook');
result = fetch(db, ['select shot, time, time_until_disrupt, v_loop, ' ...
  'ip, dipprog_dt from disruption_warning where ' ...
  'shot > 1000000 order by shot, time']);
shot = int32(cell2mat(result(:,1)));
time = cell2mat(result(:,2));
time_until_disrupt = cell2mat(result(:,3));
v_loop = cell2mat(result(:,4));
ip = cell2mat(result(:,5));
dipprog_dt = cell2mat(result(:,6));
close(db);
get_flattop_indices;
get_ip_sign;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
plot(-smear(time_until_disrupt(indices_flattop))*1.e3, ...            
v_loop(indices_flattop).*ipsign(indices_flattop),'.');
xlim([-20,0]);
ylim([0,50]);
set(gca, 'fontsize', 12);
xlabel('Time before disrupt [ms]', 'FontSize', 14);
ylabel('Loop voltage [V]', 'FontSize', 14);
title({'Loop voltage vs time before disrupt', ...
       'on Alcator C-Mod (flattop ONLY)'}, 'FontSize', 16);

set(gcf, 'PaperPositionMode', 'auto');
print -dpng vloop_plot_for_IAEA_synopsis.png;
