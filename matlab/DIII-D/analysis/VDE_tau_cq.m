cd('/fusion/projects/disruption_warning/TCV-data');

%{
fid = fopen('VDE_shotlist.txt');
VDE_shotlist = fscanf(fid, '%i');
fclose(fid);

VDE_shotlist = int32(VDE_shotlist);
nshots = length(VDE_shotlist);
taucq_80_20 = NaN(1, nshots);
z_extrema = NaN(1, nshots);

for ishot = 1:nshots;
  shot = VDE_shotlist(ishot);
  fprintf(1, 'Opening shot %i (%3i/%3i, %5.1f%%)\n', shot, ishot, nshots, ...
    100*ishot/nshots);
  [~, status] = mdsopen('magnetics', shot);
  if (~status);
    fprintf('Unable to open magnetics shot %i\n', shot);
    continue;
  end;
  [ip, status] = mdsvalue('\ip');
  if (~status);
    fprintf('Error reading Ip\n');
    continue;
  end;
  iptime = mdsvalue('dim_of(\ip)');
  zc = mdsvalue('\top.disruptions:zcentroid');
  if (~status);
    fprintf('Error reading Zcentroid\n');
    continue;
  end;
  zctime = mdsvalue('dim_of(\top.disruptions:zcentroid)');
  mdsclose;

% The timebase of the zcentroid signal is always 20 ms long, and starts 10
% ms before the disruption.  (I wrote the program that calculates the
% zcentroid signal for shots that disrupted, so that's how I know the
% timebase information.)  The plasma current on VDE's often starts to fall
% about 10 ms before the current quench (CQ), so I will determine the
% pre-VDE current by taking a 1 ms average at 20 ms before the CQ.

  tdis = mean(zctime);
  [~, indx1] = min(abs(iptime - (tdis - 20e-3 - 0.5e-3)));
  [~, indx2] = min(abs(iptime - (tdis - 20e-3 + 0.5e-3)));
  Ip0 = mean(ip(indx1:indx2));

% Find the time when Ip = 80% of Ip0.  Use linear interpolation between the
% Ip sampling times to improve resolution.

  indx = find(sign(ip(2:end) - 0.80*Ip0) ~= sign(ip(1:end-1) - 0.80*Ip0));
  indx = max(indx);
  t80 = iptime(indx) + (0.80*Ip0 - ip(indx)) * ...
    (iptime(indx+1) - iptime(indx)) / (ip(indx+1) - ip(indx));

% Find the time when Ip = 20% of Ip0.  Use linear interpolation between the
% Ip sampling times to improve resolution.

  indx = find(sign(ip(2:end) - 0.20*Ip0) ~= sign(ip(1:end-1) - 0.20*Ip0));
  indx = max(indx);
  t20 = iptime(indx) + (0.20*Ip0 - ip(indx)) * ...
    (iptime(indx+1) - iptime(indx)) / (ip(indx+1) - ip(indx));

  taucq_80_20(ishot) = t20 - t80;
  fprintf(1, '  taucq_80_20 = %5.2f ms\n', taucq_80_20(ishot)*1e3);

% Find the extrema of Zcentroid during the current quench

  indices = find(zctime >= tdis-0.005 & zctime <= tdis+0.002);
  z_extrema(ishot) = max(abs(zc(indices)));
  zmean = mean(zc(indices));
  if (zmean < 0);
    z_extrema(ishot) = -z_extrema(ishot);
  end;
  fprintf(1, '  z_extrema = %5.1f cm\n', z_extrema(ishot)*1e2);
  
end;

[~,~] = unix('chmod u+w taucq_80_20.mat');
save('taucq_80_20.mat', 'VDE_shotlist', 'taucq_80_20', 'z_extrema');
[~,~] = unix('chmod a=r taucq_80_20.mat');
%}

load('taucq_80_20.mat');

% Plot histogram of taucq_80_20
F1 = figure;
set(gca, 'fontsize', 12);
xlim([0,5]);
xlabel('\tau_{80-20}/0.6  [ms]', 'fontsize', 14);
ylabel('Probability distribution function', 'fontsize', 14);
title('\tau_{80-20}/0.6   for C-Mod VDE''s', 'fontsize', 15);
hold on;
histogram(taucq_80_20*1e3/0.6, [0 : 0.1 : 5], 'normalization', 'pdf');
hold off;
drawnow;

% Plot histogram of taucq_80_20/area
S = 0.22; % [m2] typical plasma cross-section area during flattop
F2 = figure;
set(gca, 'fontsize', 12);
xlim([0,23]);
ylim([0,0.27]);
xlabel('[\tau_{80-20}/0.6] / S  [ms/m^2]', 'fontsize', 14);
ylabel('Probability distribution function', 'fontsize', 14);
title('[\tau_{80-20}/0.6] / S   for C-Mod VDE''s', 'fontsize', 15);
hold on;
histogram(taucq_80_20*1e3/0.6/S, [0 : 0.1/0.22 : 5/0.22], ...
  'normalization', 'pdf');
hold off;
drawnow;

% Plot z_extrema vs shot #.  This shows that VDE's can go up or down, and
% doesn't seem to change over the years.

F3 = figure;
set(gca, 'fontsize', 12);
xlim([0,length(VDE_shotlist)]);
ylim([-0.4, 0.4]);
xlabel('VDE shotlist index #', 'fontsize', 14);
ylabel('Zcentroid extrema [m]', 'fontsize', 14);
title('Zcentroid extrema  for C-Mod VDE''s', 'fontsize', 15);
hold on;
plot(z_extrema, '.b');
hold off;
drawnow;

% Plot [taucq_80_20 / S] vs z_extrema.

F4 = figure;
set(gca, 'fontsize', 12);
xlim([-0.4, 0.4]);
ylim([0, 23e-3]);
xlabel('Zcentroid extrema [m]', 'fontsize', 14);
ylabel('[\tau_{80-20}/0.6] / S  [ms/m^2]', 'fontsize', 14);
title('[\tau_{80-20}/0.6]/S vs Z extrema  for C-Mod VDE''s', 'fontsize', 15);
hold on;
plot(z_extrema, taucq_80_20/0.6/S, '.b');
hold off;
drawnow;

data = readmatrix('cmod_ITPA_CQ_Data.txt');
ITPA_cmod_taucq_80_20_over_06_over_area = data(:,2);
clearvars data;

F5 = figure;
xlim([0,23]);
ylim([0,0.27]);
xlabel('[\tau_{80-20}/0.6] / S  [ms/m^2]', 'fontsize', 14);
ylabel('Probability distribution function', 'fontsize', 14);
title('[\tau_{80-20}/0.6] / S   for C-Mod MD''s and VDE''s', 'fontsize', 15);
hold on;

histogram(taucq_80_20*1e3/0.6/S, [0 : 0.1/0.22 : 5/0.22], ...
  'normalization', 'pdf', 'facecolor', [0.000, 0.447, 0.741]);

histogram(ITPA_cmod_taucq_80_20_over_06_over_area, [0: 0.1/0.22 : 50/0.22], ...
  'normalization', 'pdf', 'facecolor', 'red', ...
  'facealpha', 0.6);

legend('C-Mod VDE''s only', ...
  [newline 'All C-Mod disruptions' newline '   in ITPA database']);
hold off;
drawnow;
