% retrieve_all_params;

get_ip_sign;
define_indices;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none', ...
  'PaperPositionMode', 'auto');

plot([-0.030,0], [0,0], 'k', 'linewidth', 1.5);
xlim([-0.030,0]); ylim([-500,100]);
set(gca, 'linewidth', 1.5, 'fontsize', 12)
hold on
plot(-time_until_disrupt(indices_flattop_disrupt_in_flattop),...
  ip_error(indices_flattop_disrupt_in_flattop)/1.e3 .* ...
  ipsign(indices_flattop_disrupt_in_flattop), '.');

xlabel('Time before disrupt [s]', 'FontSize', 14);
ylabel('Ip error [kA]', 'FontSize', 14);
title('C-Mod: I_p error', 'FontSize', 16);
%text(0.65,0.85,{'all flattop disruptions'}, 'fontsize', 16, ... 
%  'units', 'normalized', 'edgecolor', 'k', 'margin', 5);
print('cmod_ip_error_vs_time_disrupt_unsmeared.png', '-dpng');

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none', ...
  'PaperPositionMode', 'auto');

plot([-0.030,0], [0,0], 'k', 'linewidth', 1.5);
xlim([-0.030,0]); ylim([-500,100]);
set(gca, 'linewidth', 1.5, 'fontsize', 12)
hold on
plot(-smear(time_until_disrupt(indices_flattop_disrupt_in_flattop),.0005),...
  ip_error(indices_flattop_disrupt_in_flattop)/1.e3 .* ...
  ipsign(indices_flattop_disrupt_in_flattop), '.');

xlabel('Time before disrupt [s]', 'FontSize', 14);
ylabel('Ip error [kA]', 'FontSize', 14);
title('C-Mod: I_p error', 'FontSize', 16);
%text(0.65,0.85,{'all flattop disruptions'}, 'fontsize', 16, ... 
%  'units', 'normalized', 'edgecolor', 'k', 'margin', 5);
print('cmod_ip_error_vs_time_disrupt.png', '-dpng');
