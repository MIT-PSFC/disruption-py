% Get data
shot = 175697;
efit_tree  = select_efit_trees(shot, '', 'DIS');
timebase = [1:0.005:3.1];
disp("Getting Kevin's peaking factors ...")
[Te_pf, ne_pf, Rad_CVA, Rad_XDIV] = get_peaking_factors_d3d(shot,timebase);
disp("Getting Ale's peaking factors ...")
Profiles = PF_plot(shot,[min(timebase),max(timebase)],'custom',efit_tree{1})

% Plot
figure();
subplot(4,1,1); hold on;
plot(timebase, Te_pf,'DisplayName','Kevin')
plot(Profiles{1}{2}.time,Profiles{1}{1},'DisplayName','Ale')
ylabel('Te PF')
lgd = legend();
lgd.('Interpreter') = 'latex';
subplot(4,1,2); hold on;
plot(timebase, ne_pf,'DisplayName','Kevin')
plot(Profiles{2}{2}.time,Profiles{2}{1},'DisplayName','Ale')
ylabel('ne PF')
lgd = legend();
lgd.('Interpreter') = 'latex';
subplot(4,1,3); hold on;
plot(timebase, Rad_CVA,'DisplayName','Kevin')
plot(Profiles{3}{2}.time,Profiles{3}{1},'DisplayName','Ale')
ylabel('Rad\_CVA')
lgd = legend();
lgd.('Interpreter') = 'latex';
subplot(4,1,4); hold on;
plot(timebase, Rad_XDIV,'DisplayName','Kevin')
plot(Profiles{4}{2}.time,Profiles{4}{1},'DisplayName','Ale')
ylabel('Rad\_XDIV')
lgd = legend();
lgd.('Interpreter') = 'latex';
xlabel('Time [s]')
title(['Shot #' num2str(shot)])

