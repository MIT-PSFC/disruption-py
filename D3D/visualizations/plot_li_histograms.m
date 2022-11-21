dbconn = set_database('d3drdb');

if exist('shot','var') == 0;
  result = fetch(dbconn,['select shot from ' ...
    'disruption_warning order by shot, time']);
  shot = int32(cell2mat(result));
end;
if exist('time','var') == 0;
  result = fetch(dbconn,['select time from ' ...
    'disruption_warning order by shot, time']);
  time = cell2mat(result);
end;
if exist('time_until_disrupt','var') == 0;
  result = fetch(dbconn,['select time_until_disrupt from ' ...
    'disruption_warning order by shot, time']);
  time_until_disrupt = cell2mat(result);
end;
if exist('ip','var') == 0;
  result = fetch(dbconn,['select ip from ' ...
    'disruption_warning order by shot, time']);
  ip = cell2mat(result);
end;
if exist('ip_error','var') == 0;
  result = fetch(dbconn,['select ip_error from ' ...
    'disruption_warning order by shot, time']);
  ip_error = cell2mat(result);
end;
if exist('li','var')==0;
  result = fetch(dbconn,['select li from ' ...
    'disruption_warning order by shot, time']);
  li = cell2mat(result);
end;
if exist('dipprog_dt','var')==0;
  result = fetch(dbconn,['select dipprog_dt from ' ...
    'disruption_warning order by shot, time']);
  dipprog_dt = cell2mat(result);
end;
if exist('power_supply_railed','var')==0;
  result = fetch(dbconn,['select power_supply_railed from ' ...
    'disruption_warning order by shot, time']);
  power_supply_railed = cell2mat(result);
end;
if exist('intentional_disruption','var')==0;
  result = fetch(dbconn,['select intentional_disruption from ' ...
    'disruption_warning order by shot, time']);
  intentional_disruption = cell2mat(result);
end;
if exist('other_hardware_failure','var')==0;
  result = fetch(dbconn,['select other_hardware_failure from ' ...
    'disruption_warning order by shot, time']);
  other_hardware_failure = cell2mat(result);
end;
close(dbconn);

define_indices;
ii_far = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 1.1));
ii_10 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 1.0 & time_until_disrupt <= 1.1));
ii_09 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.9 & time_until_disrupt <= 1.0));
ii_08 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.8 & time_until_disrupt <= 0.9));
ii_07 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.7 & time_until_disrupt <= 0.8));
ii_06 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.6 & time_until_disrupt <= 0.7));
ii_05 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.5 & time_until_disrupt <= 0.6));
ii_04 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.4 & time_until_disrupt <= 0.5));
ii_03 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.3 & time_until_disrupt <= 0.4));
ii_02 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.2 & time_until_disrupt <= 0.3));
ii_01 = intersect(indices_flattop_disrupt_in_flattop, ...
  find(time_until_disrupt >= 0.1 & time_until_disrupt <= 0.2));

xaxis = [-1.05, -0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15];

bins = 0 : .02 : 2;
bin_centers = bins(1:end-1) + diff(bins)/2;

h0 = histogram(li(ii_far), bins, 'normalization', 'probability' , ...
  'facecolor', 'green');
set(gca, 'fontsize', 12);
xlabel('l_i', 'fontsize', 14);
ylabel('probability', 'fontsize', 14);
title('Histograms of l_i at various times before disrupt', 'fontsize', 16);
hold on;

h10 = histogram(li(ii_10), bins, 'normalization', 'probability' , ...
  'facecolor', 'black');
h9 = histogram(li(ii_09), bins, 'normalization', 'probability' , ...
  'facecolor', 'blue');
h8 = histogram(li(ii_08), bins, 'normalization', 'probability' , ...
  'facecolor', 'red');
h7 = histogram(li(ii_07), bins, 'normalization', 'probability' , ...
  'facecolor', 'yellow');
h6 = histogram(li(ii_06), bins, 'normalization', 'probability' , ...
  'facecolor', 'green');
h5 = histogram(li(ii_05), bins, 'normalization', 'probability' , ...
  'facecolor', 'blue');
h4 = histogram(li(ii_04), bins, 'normalization', 'probability' , ...
  'facecolor', 'red');
h3 = histogram(li(ii_03), bins, 'normalization', 'probability' , ...
  'facecolor', 'yellow');
h2 = histogram(li(ii_02), bins, 'normalization', 'probability' , ...
  'facecolor', 'green');
h1 = histogram(li(ii_01), bins, 'normalization', 'probability' , ...
  'facecolor', 'blue');

% Re-plot li(ii_far) so it can be seen on top of all the others.

histogram(li(ii_far), bins, 'normalization', 'probability' , ...
  'facecolor', 'green');

hold off;

mean_li = NaN(1,11);

mean_li(1) = xcentroid(bin_centers, h0.Values);
mean_li(2) = xcentroid(bin_centers, h10.Values);
mean_li(3) = xcentroid(bin_centers, h9.Values);
mean_li(4) = xcentroid(bin_centers, h8.Values);
mean_li(5) = xcentroid(bin_centers, h7.Values);
mean_li(6) = xcentroid(bin_centers, h6.Values);
mean_li(7) = xcentroid(bin_centers, h5.Values);
mean_li(8) = xcentroid(bin_centers, h4.Values);
mean_li(9) = xcentroid(bin_centers, h3.Values);
mean_li(10) = xcentroid(bin_centers, h2.Values);
mean_li(11) = xcentroid(bin_centers, h1.Values);

figure;
plot(xaxis, mean_li(2:end), 'sb'); xlim([-1.5,0]);
set(gca, 'fontsize', 12);
xlabel('Time until disrupt [s]', 'fontsize', 14);
ylabel('mean(l_i)', 'fontsize', 14);
title('Mean of l_i distribution vs time before disrupt', ...
  'fontsize', 16);
hold on;
plot([-1.5,-1.1], [mean_li(1), mean_li(1)], '-b');

difffrac = NaN(1,10);

indices = find(h10.Values > h0.Values & bin_centers > 1.0);
difffrac(1) = sum(h10.Values(indices) - h0.Values(indices));
indices = find(h9.Values > h0.Values & bin_centers > 1.0);
difffrac(2) = sum(h9.Values(indices) - h0.Values(indices));
indices = find(h8.Values > h0.Values & bin_centers > 1.0);
difffrac(3) = sum(h8.Values(indices) - h0.Values(indices));
indices = find(h7.Values > h0.Values & bin_centers > 1.0);
difffrac(4) = sum(h7.Values(indices) - h0.Values(indices));
indices = find(h6.Values > h0.Values & bin_centers > 1.0);
difffrac(5) = sum(h6.Values(indices) - h0.Values(indices));
indices = find(h5.Values > h0.Values & bin_centers > 1.0);
difffrac(6) = sum(h5.Values(indices) - h0.Values(indices));
indices = find(h4.Values > h0.Values & bin_centers > 1.0);
difffrac(7) = sum(h4.Values(indices) - h0.Values(indices));
indices = find(h3.Values > h0.Values & bin_centers > 1.0);
difffrac(8) = sum(h3.Values(indices) - h0.Values(indices));
indices = find(h2.Values > h0.Values & bin_centers > 1.0);
difffrac(9) = sum(h2.Values(indices) - h0.Values(indices));
indices = find(h1.Values > h0.Values & bin_centers > 1.0);
difffrac(10) = sum(h1.Values(indices) - h0.Values(indices));

figure;
plot(xaxis, difffrac, 'sb'); xlim([-1.10,0]);
set(gca, 'fontsize', 12);
xlabel('Time until disrupt [s]', 'fontsize', 14);
ylabel('fractional shift', 'fontsize', 14);
title('Fractional shift of l_i distribution vs time before disrupt', ...
  'fontsize', 16);
