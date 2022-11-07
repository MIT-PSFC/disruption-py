if (exist('time','var')==0 || exist('time_until_disrupt','var')==0 || ...
    exist('v_loop','var')==0);
  fprintf(1,['Call the routine ''retrieve_all_params'' before calling ' ...
             'this routine.\n']);
  return;
end;

get_ip_sign;
get_flattop_indices;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
plot(-time_until_disrupt(indices_flattop), ...         
  v_loop(indices_flattop).*ipsign(indices_flattop),'.');
xlim([-1,0]);
ylim([0,9]);
xlabel('Time before disrupt [s]', 'FontSize', 14);
ylabel('Loop voltage (processed for sign of Ip) [V]', 'FontSize', 14);
title('Loop voltage vs time before disrupt (flattop ONLY)', 'FontSize', 16);
 
set(gcf, 'PaperPositionMode', 'auto');
print -dpng vloop_plot_1.png;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
plot(-smear(time_until_disrupt(indices_flattop))*1.e3, ...         
  v_loop(indices_flattop).*ipsign(indices_flattop),'.');
xlim([-50,0]);
ylim([0,9]);
xlabel('Time before disrupt [ms]', 'FontSize', 14);
ylabel('Loop voltage (processed for sign of Ip) [V]', 'FontSize', 14);
title('Loop voltage vs time before disrupt (flattop ONLY)', 'FontSize', 16);
 
set(gcf, 'PaperPositionMode', 'auto');
print -dpng vloop_plot_2.png;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')

indices_nodisrupt = find(isnan(time_until_disrupt));
indices_nodisrupt_flattop = intersect(indices_nodisrupt, indices_flattop);

plot(time(indices_nodisrupt_flattop), ...         
  v_loop(indices_nodisrupt_flattop).*ipsign(indices_nodisrupt_flattop),'.');
xlim([0,2]);
ylim([0,9]);
xlabel('Time [s]', 'FontSize', 14);
ylabel('Loop voltage (processed for sign of Ip) [V]', 'FontSize', 14);
title('Loop voltage for non-disruptive shots (flattop ONLY)', 'FontSize', 16);
 
set(gcf, 'PaperPositionMode', 'auto');
print -dpng vloop_plot_3.png;
