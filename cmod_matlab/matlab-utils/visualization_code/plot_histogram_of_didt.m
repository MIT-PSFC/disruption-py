db = set_database('logbook');
result = fetch(db,['select shot,time,ip,dipprog_dt from ' ...
                   'disruption_warning order by shot,time']);
close(db);

shots = int32(cell2mat(result(:,1)));
times = cell2mat(result(:,2));
ip = cell2mat(result(:,3));
didt = cell2mat(result(:,4));

% Take account of the Ip polarity.  First I must create an array containing
% the Ip polarity for each value of didt.

polarity = zeros(length(didt),1);
unique_shots = unique(shots);
for i = 1:length(unique_shots);
  indices = find(shots == unique_shots(i));
  polarity(indices) = sign(sum(ip(indices)));
end;

didt = didt .* polarity;

% Calculate and plot the histogram of dI/dt values for each time slice

[timevals, indices, jndices] = unique(times);
ntimes = length(timevals);

didt_bin_width = 1.e6 * 0.2; % 0.2 MA/s (or kA/ms)
didt_bins = [-1.e7 : didt_bin_width : +1.e7];
didt_bin_edges = [didt_bins - 0.5 * didt_bin_width, +Inf];
nbins = length(didt_bins);

hist2d = zeros(ntimes, nbins);
for itime = 1:ntimes;
  didtvals = didt(find(abs(times - timevals(itime)) < 0.010));
  [hist1d, dummy] = histcounts(didtvals, didt_bin_edges, 'Normalization', ...
                               'probability');
  hist2d(itime, :) = hist1d;
end;

figure(2);
for itime = 1:ntimes;
  bar(didt_bins/1.e6, hist2d(itime, :), 1);
  xlim([-6, +6]);
  ylim([0,1]);
  xlabel('programmed dIp/dt [MA/s or kA/ms]');
  ylabel('probability');
  text(0.05, 0.90, ['Time = ' num2str(timevals(itime),'%5.3f') ' s'], ...
                    'units', 'normalized');
  fprintf(1,'Time = %5.3f s\n', timevals(itime));
  pause;
end;
