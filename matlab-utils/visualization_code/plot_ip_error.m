if (exist('time','var')==0 || exist('time_until_disrupt','var')==0 || ...
    exist('ip','var')==0 || exist('ip_error','var')==0);
  fprintf(1,['Call the routine ''retrieve_all_params'' before calling ' ...
             'this routine.\n']);
  return;
end;

get_ip_sign;
get_flattop_indices;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
plot(-time_until_disrupt(indices_flattop), ...         
  ip_error(indices_flattop)/1.e3.*ipsign(indices_flattop),'.');
xlim([-1,0]);
ylim([-700,200]);
xlabel('Time before disrupt [s]', 'FontSize', 14);
ylabel('Ip error (processed for sign of Ip) [kA]', 'FontSize', 14);
title('Ip error vs time before disrupt (flattop ONLY)', 'FontSize', 16);
 
set(gcf, 'PaperPositionMode', 'auto');
%print -dpng iperror_plot_1.png;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
plot(-smear(time_until_disrupt(indices_flattop))*1.e3, ...         
  ip_error(indices_flattop)/1.e3.*ipsign(indices_flattop),'.');
xlim([-50,0]);
ylim([-700,200]);
xlabel('Time before disrupt [ms]', 'FontSize', 14);
ylabel('Ip error (processed for sign of Ip) [kA]', 'FontSize', 14);
title('Ip error vs time before disrupt (flattop ONLY)', 'FontSize', 16);
 
set(gcf, 'PaperPositionMode', 'auto');
%print -dpng iperror_plot_2.png;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')

indices_nodisrupt = find(isnan(time_until_disrupt));
indices_nodisrupt_flattop = intersect(indices_nodisrupt, indices_flattop);

plot(time(indices_nodisrupt_flattop), ...         
  ip_error(indices_nodisrupt_flattop)/1.e3 .* ...
  ipsign(indices_nodisrupt_flattop),'.');
xlim([0,2]);
ylim([-700,200]);
xlabel('Time [s]', 'FontSize', 14);
ylabel('Ip error (processed for sign of Ip) [kA]', 'FontSize', 14);
title('Ip error for non-disruptive shots (flattop ONLY)', 'FontSize', 16);
 
set(gcf, 'PaperPositionMode', 'auto');
%print -dpng iperror_plot_3.png;
