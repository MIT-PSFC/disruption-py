% retrieve_all_params;

get_ip_sign;
define_indices;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none', ...
  'PaperPositionMode', 'auto');

plot([0,2], [0,0], 'k', 'linewidth', 1.5);
xlim([0,2]); ylim([-500,100]);
set(gca, 'linewidth', 1.5, 'fontsize', 12)
hold on
plot(smear(time(indices_flattop_no_disrupt),0.000),...
  ip_error(indices_flattop_no_disrupt)/1.e3 .* ...
  ipsign(indices_flattop_no_disrupt), '.');

xlabel('Time [s]', 'FontSize', 14);
ylabel('Ip error [kA]', 'FontSize', 14);
title('C-Mod: I_p error', 'FontSize', 16);
%text(0.65,0.85,{'all flattop disruptions'}, 'fontsize', 16, ... 
%  'units', 'normalized', 'edgecolor', 'k', 'margin', 5);
print('cmod_ip_error_vs_time_no_disrupt.png', '-dpng');
