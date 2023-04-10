if (exist('time','var')==0 || exist('time_until_disrupt','var')==0 || ...
    exist('ip','var')==0 || exist('ip_error','var')==0);
  fprintf(1,['Call the routine ''retrieve_all_params_d3d'' before calling ' ...
             'this routine.\n']);
  return;
end;

get_ip_sign;
get_flattop_indices;

% Get indices of time of disruption for all shots that disrupted during
% flattop:
qq = intersect(find(time_until_disrupt==0), indices_flattop);

% Therefore, shot(qq) are all shots that disrupted during flattop
% Get all the indices for shots that disrupted during flattop
qq = find(ismember(shot, shot(qq)));

% Finally, get the subset of the qq time slices that occur during flattop
qq = intersect(qq, indices_flattop);

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
plot(-time_until_disrupt(qq), ip_error(qq)/1.e3.*ipsign(qq),'.');
xlim([-1,0]);
ylim([-700,200]);
xlabel('Time before disrupt [s]', 'FontSize', 16);
ylabel('Ip error [kA]', 'FontSize', 16);
%title('Ip error vs time before disrupt (flattop ONLY)', 'FontSize', 16);
 
set(gcf, 'PaperPositionMode', 'auto');
print -dpng iperror_d3d_plot_1.png;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')
plot(-smear(time_until_disrupt(qq)), ...         
  ip_error(qq)/1.e3.*ipsign(qq),'.');
xlim([-0.25,0]);
ylim([-700,200]);
xlabel('Time before disrupt [s]', 'FontSize', 16);
ylabel('Ip error [kA]', 'FontSize', 16);
%title('Ip error vs time before disrupt (flattop ONLY)', 'FontSize', 16);
 
set(gcf, 'PaperPositionMode', 'auto');
print -dpng iperror_d3d_plot_2.png;

figure('Name', ' ', 'NumberTitle','off','Menubar','none','Toolbar','none')

indices_nodisrupt = find(isnan(time_until_disrupt));
indices_nodisrupt_flattop = intersect(indices_nodisrupt, indices_flattop);

plot(time(indices_nodisrupt_flattop), ...         
  ip_error(indices_nodisrupt_flattop)/1.e3 .* ...
  ipsign(indices_nodisrupt_flattop),'.');
xlim([0,9]);
ylim([-700,200]);
xlabel('Time [s]', 'FontSize', 16);
ylabel('Ip error [kA]', 'FontSize', 16);
%title('Ip error for non-disruptive shots (flattop ONLY)', 'FontSize', 16);
 
set(gcf, 'PaperPositionMode', 'auto');
print -dpng iperror_d3d_plot_3.png;
