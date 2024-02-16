%hist3([-time_until_disrupt(indices_flattop_disrupt_in_flattop),beta_n(indices_flattop_disrupt_in_flattop)])
index = find(time_until_disrupt(indices_flattop_disrupt_in_flattop)<0.3 & beta_n(indices_flattop_disrupt_in_flattop)>=0 & beta_n(indices_flattop_disrupt_in_flattop)<=4);
tud = time_until_disrupt(indices_flattop_disrupt_in_flattop);
bn = beta_n(indices_flattop_disrupt_in_flattop);
x = tud(index);
y = bn(index);
ctrs = {0.025:0.025:0.275 0.25:.25:3.75}
hist3([x,y],'Ctrs',ctrs)
xlabel('time before disrupt [s]');
ylabel('\beta_n');
%zlim([0,6500]);
hold on;
N = hist3([x,y],'Ctrs',ctrs);
N_pcolor = N';
N_pcolor(size(N_pcolor,1)+1,size(N_pcolor,2)+1) = 0;
xl = linspace(min(x),max(x),size(N_pcolor,2)); % Columns of N_pcolor
yl = linspace(min(y),max(y),size(N_pcolor,1)); % Rows of N_pcolor
h = pcolor(xl,yl,N_pcolor);
colormap('hot') % Change color scheme 
colorbar % Display colorbar
h.ZData = -max(N_pcolor(:))*ones(size(N_pcolor));
ax = gca;
ax.ZTick(ax.ZTick < 0) = [];
title({'\beta_n Histogram and Intensity Map','flattop disruptions time slices'});

%betan
X = [x,y];
ctrs = {0.025:0.025:0.275 0.25:.25:3.75}
hist3(X,'CDataMode','auto','FaceColor','interp','Ctrs',ctrs)
colorbar
view(2);
xlabel('time before disrupt [s]')
ylabel('\beta_n')
title('DIII-D flattop disruptions time slices')

%qstar
index = find(time_until_disrupt(indices_flattop_disrupt_in_flattop)>0.005 & time_until_disrupt(indices_flattop_disrupt_in_flattop)<0.3 & qstar(indices_flattop_disrupt_in_flattop)<12);
tud = time_until_disrupt(indices_flattop_disrupt_in_flattop);
qs = qstar(indices_flattop_disrupt_in_flattop);

x = tud(index);
y = qs(index);

X = [x,y];
ctrs = {0.025:0.025:0.275 0.25:.4:12}
hist3(X,'CDataMode','auto','FaceColor','interp','Ctrs',ctrs)
colorbar
view(2);
xlabel('time before disrupt [s]')
ylabel('q^*')
title('DIII-D flattop disruptions time slices')

%kappa_area
index = find(time_until_disrupt(indices_flattop_disrupt_in_flattop)>0.005 & time_until_disrupt(indices_flattop_disrupt_in_flattop)<0.3);
tud = time_until_disrupt(indices_flattop_disrupt_in_flattop);
ka = kappa_area(indices_flattop_disrupt_in_flattop);

x = tud(index);
y = ka(index);

X = [x,y];
ctrs = {0.025:0.025:0.275 1:.05:2}
hist3(X,'CDataMode','auto','FaceColor','interp','Ctrs',ctrs)
colorbar
view(2);
xlabel('time before disrupt [s]')
ylabel('\kappa_{a}')
title('DIII-D flattop disruptions time slices')

%triangularity
index = find(time_until_disrupt(indices_flattop_disrupt_in_flattop)>0.005 & time_until_disrupt(indices_flattop_disrupt_in_flattop)<0.3 & delta(indices_flattop_disrupt_in_flattop) > 0);
tud = time_until_disrupt(indices_flattop_disrupt_in_flattop);
d = delta(indices_flattop_disrupt_in_flattop);

x = tud(index);
y = d(index);

X = [x,y];
ctrs = {0.025:0.025:0.275 0:.05:1};
hist3(X,'CDataMode','auto','FaceColor','interp','Ctrs',ctrs)
colorbar
view(2);
xlabel('time before disrupt [s]')
ylabel('\delta')
title('DIII-D flattop disruptions time slices')

%betap
index = find(time_until_disrupt(indices_flattop_disrupt_in_flattop)>0.005 & time_until_disrupt(indices_flattop_disrupt_in_flattop)<0.3 & beta_p(indices_flattop_disrupt_in_flattop) > 0);
tud = time_until_disrupt(indices_flattop_disrupt_in_flattop);
bp = beta_p(indices_flattop_disrupt_in_flattop);

x = tud(index);
y = bp(index);

X = [x,y];
ctrs = {0.025:0.025:0.275 0.06:0.06:2};
hist3(X,'CDataMode','auto','FaceColor','interp','Ctrs',ctrs)
colorbar
view(2);
xlabel('time before disrupt [s]')
ylabel('\beta_p')
title('DIII-D flattop disruptions time slices')

%q95
index = find(time_until_disrupt(indices_flattop_disrupt_in_flattop)>0.005 & time_until_disrupt(indices_flattop_disrupt_in_flattop)<0.3 & q95(indices_flattop_disrupt_in_flattop) < 10);
tud = time_until_disrupt(indices_flattop_disrupt_in_flattop);
q = q95(indices_flattop_disrupt_in_flattop);

x = tud(index);
y = q(index);

X = [x,y];
ctrs = {0.025:0.025:0.275 1.9:0.12:7};
hist3(X,'CDataMode','auto','FaceColor','interp','Ctrs',ctrs)
colorbar
view(2);
xlabel('time before disrupt [s]')
ylabel('q95')
title('DIII-D flattop disruptions time slices')

%Greenwald
index = find(time_until_disrupt(indices_flattop_disrupt_in_flattop)>0.005 & time_until_disrupt(indices_flattop_disrupt_in_flattop)<0.3 & Greenwald_fraction(indices_flattop_disrupt_in_flattop) < 2);
tud = time_until_disrupt(indices_flattop_disrupt_in_flattop);
gf = Greenwald_fraction(indices_flattop_disrupt_in_flattop);

x = tud(index);
y = gf(index);

X = [x,y];
ctrs = {0.025:0.025:0.275 0:0.05:1};
hist3(X,'CDataMode','auto','FaceColor','interp','Ctrs',ctrs)
colorbar
view(2);
xlabel('time before disrupt [s]')
ylabel('ne/nG')
title('DIII-D flattop disruptions time slices')
%print('/fusion/projects/disruption_warning/matlab_programs/SPARC\ hist/Green.png', '-dpng');

