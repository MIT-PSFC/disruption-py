function [shotlist, probability] = plot_stats(averaging_window_size, ...
          shot_range);

if (exist('averaging_window_size') && ...
    ~isempty(averaging_window_size) && ...
    (averaging_window_size >= 1));
  window = averaging_window_size;
else;
  window = 10;
end;

if (exist('shot_range') && ...
    ~isempty(shot_range) && ...
    numel(shot_range) == 2);
  shotrange = shot_range;
else;
  shotrange = [100000, 161650];
end;

shotlist = zeros(1,161774);
probability = zeros(1,161774);

db = set_database('d3drdb');
ps = db.prepareStatement('select shot from disruptions order by shot');
rs = ps.executeQuery();
while (rs.next());
  disruption_shot = rs.getObject(1);
  shotlist(disruption_shot) = 1;
end;
db.close();
ps.close();
rs.close();

if (mod(window,2)==1);
  pre_window = (window-1)/2;
else;
  pre_window = window/2;
end;
post_window = window - pre_window - 1;

for i = shotrange(1):shotrange(2);
  probability(i) = mean(shotlist(max([i- pre_window,  1]) : ...
                                 min([i+post_window,end])));
end;

plot([shotrange(1):shotrange(2)], probability(shotrange(1):shotrange(2)), ...
     '-');
set(gcf, 'menubar', 'none');
xlabel('DIII-D shot number', 'fontsize', 12);
ylabel(['Disruption probability (averaged over ' num2str(window) ...
        ' shots)'], 'fontsize', 12);
title('Disruption probability history', 'fontsize', 14);
ylim([0,1]);
ytickvals = [0 : 0.25 : 1];
yticklabels = cell(1,numel(ytickvals));
for itick = 1 : numel(ytickvals);
  yticklabels(itick) = cellstr(sprintf('%4.2f',ytickvals(itick)));
end;
set(gca, 'YTick', ytickvals, 'YTickLabel', yticklabels, ...
         'YMinorTick', 'on');

xtickvals = get(gca, 'XTick');
xticklabels = cell(1,numel(xtickvals));
for itick = 1 : numel(xtickvals);
  xticklabels(itick) = cellstr(sprintf('%6i',xtickvals(itick)));
end;
set(gca, 'XTickLabel', xticklabels);

if numel(xtickvals > 5);
  xtickvals = input('Enter new set of X tick values: ');
  if ~isempty(xtickvals);
    if (numel(xtickvals) <= 5);
      xminorticks = 'on';
    else;
      xminorticks = 'off';
    end;
    xticklabels = cell(1,numel(xtickvals));
    for itick = 1 : numel(xtickvals);
      xticklabels(itick) = cellstr(sprintf('%6i',xtickvals(itick)));
    end;
    set(gca, 'XTick', xtickvals, 'XTickLabel', xticklabels, ...
             'XMinorTick', xminorticks);
  end;
end;
