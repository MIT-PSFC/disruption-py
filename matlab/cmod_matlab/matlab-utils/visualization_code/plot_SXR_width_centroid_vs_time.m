fprintf(1, ['Before running this program, all the relevant data \n' ...
            'must be read out of the disruption warning database.\n']);

define_indices;
ii = indices_flattop_no_disrupt;
jj = indices_flattop_disrupt_in_flattop;

times = [0.01 : 0.01 : 10.00];
ntimes = length(times);
centroid = NaN(size(times));

for it = 1:ntimes;
  kk = intersect(jj, find(time_until_disrupt > times(it) - 0.005 & ...
                          time_until_disrupt < times(it) + 0.005));
  centroid(it) = mean(SXR_peaking(kk), 'omitnan');
end;

figure;
plot(-times, centroid, '.-b');
xlim([-5,0]);
ylim([1,1.5]);
set(gca, 'fontsize', 12);
xlabel('Time until disrupt [s]', 'fontsize', 14);
ylabel('<SXR peaking>', 'fontsize', 14);
title('<SXR peaking> vs time until disrupt', 'fontsize', 15);

times = [1.0 : 0.1 : 9.0];
ntimes = length(times);
centroid = NaN(size(times));

for it = 1:ntimes;
  kk = intersect(ii, find(time > times(it) - 0.050 & ...
                          time < times(it) + 0.050));
  centroid(it) = mean(SXR_peaking(kk), 'omitnan');
end;

figure;
plot(times, centroid, '.-b');
xlim([0,10]);
ylim([1,1.5]);
set(gca, 'fontsize', 12);
xlabel('Time [s]', 'fontsize', 14);
ylabel('<SXR peaking>', 'fontsize', 14);
title('<SXR peaking>, non-disruptions', 'fontsize', 15);
