if (~exist('already_have_data'));
  db = set_database('d3drdb'); % My routine to connect to database

  ps = db.prepareStatement('select count(1) as "nshots" from disruptions');
  rs = ps.executeQuery();
  rs.next();
  ndisrupt = rs.getObject(1);

  shots_disrupt = zeros(1, ndisrupt);
  ip_disrupt = zeros(1, ndisrupt);
  betan_disrupt = zeros(1, ndisrupt);
  t_disrupt = zeros(1, ndisrupt);

  ps = db.prepareStatement(['select shot,ip,betan,t_disrupt ' ...
    'from disruptions order by shot']);
  rs = ps.executeQuery();

  for ishot = 1:ndisrupt;
    rs.next();
    shots_disrupt(ishot) = rs.getObject(1);
    ip_disrupt(ishot) = rs.getObject(2);

    if isempty(rs.getObject(3));
      betan_disrupt(ishot) = NaN;
    else;
      betan_disrupt(ishot) = rs.getObject(3);
    end;

    t_disrupt(ishot) = rs.getObject(4);
  end;

  ps = db.prepareStatement(['select count(1) as "nshots" from summaries ' ...
    'where betanmax is not null and ' ...
    'shot >= ' num2str(shots_disrupt(1),'%i') ]);
  rs = ps.executeQuery();
  rs.next();
  nshots = rs.getObject(1);

  shots_summaries = zeros(1, nshots);
  ip_summaries = zeros(1, nshots);
  betanmax = zeros(1, nshots);
  t_betanmax = zeros(1, nshots);

  ps = db.prepareStatement(['select shot,ip,betanmax,t_betanmax ' ...
    'from summaries ' ...
    'where betanmax is not null ' ...
    'and shot >= ' num2str(shots_disrupt(1),'%i') ...
    'order by shot']);
  rs = ps.executeQuery();

  for ishot = 1:nshots;
    rs.next();
    shots_summaries(ishot) = rs.getObject(1);

    if isempty(rs.getObject(2));
      ip_summaries(ishot) = NaN;
    else;
      ip_summaries(ishot) = rs.getObject(2);
    end;

    betanmax(ishot) = rs.getObject(3);
    t_betanmax(ishot) = rs.getObject(4);
  end;

  db.close();
  ps.close();
  rs.close();
  clear db ps rs ishot ans;
  already_have_data = [];
end;

[shots_disrupt_summaries, index_disrupt_summaries, ...
  index_disrupt_disruptions] = intersect(shots_summaries, shots_disrupt);
[shots_non_disrupt_summaries, index_non_disrupt_summaries] = ...
  setdiff(shots_summaries, shots_disrupt);

betanmax_disrupt = betanmax(index_disrupt_summaries);
betanmax_non_disrupt = betanmax(index_non_disrupt_summaries);

dx_bin = 0.2;
xbins = [0 : dx_bin : 6-dx_bin] + dx_bin/2;

index_valid_range = find(betanmax >= 0 & betanmax <= 6);
figure(1);
set(gcf, 'menubar', 'none');
ybins_all = hist(betanmax(index_valid_range), xbins);
bar(xbins, ybins_all, 'b');
xtickvals = get(gca, 'XTick');
xticklabels = cell(1,numel(xtickvals));
for itick = 1 : numel(xtickvals);
  xticklabels(itick) = cellstr(sprintf('%3.1f',xtickvals(itick)));
end;
set(gca, 'XTickLabel', xticklabels);
ylimits = ylim;
xlabel('\beta_Nmax', 'fontsize', 12);
ylabel('Number of shots', 'fontsize', 12);
title('Histogram of \beta_Nmax for all shots', 'fontsize', 14);

index_valid_range = find(betanmax_disrupt >= 0 & ...
  betanmax_disrupt <= 6);
figure(2);
set(gcf, 'menubar', 'none');
ybins_disrupt = hist(betanmax_disrupt(index_valid_range), xbins);
bar(xbins, ybins_disrupt, 'FaceColor', 'red');
ylim(ylimits);
xtickvals = get(gca, 'XTick');
xticklabels = cell(1,numel(xtickvals));
for itick = 1 : numel(xtickvals);
  xticklabels(itick) = cellstr(sprintf('%3.1f',xtickvals(itick)));
end;
set(gca, 'XTickLabel', xticklabels);
xlabel('\beta_Nmax', 'fontsize', 12);
ylabel('Number of shots', 'fontsize', 12);
title('Histogram of \beta_Nmax for disruptive shots', 'fontsize', 14);

figure(3);
set(gcf, 'menubar', 'none');
bar(xbins, ybins_all, 'b', 'EdgeColor', 'black');
hold on;
bar(xbins+0.05, ybins_disrupt, 'r', 'EdgeColor', 'black');
hold off;
xlim([0,6]);
ylim(ylimits);
xtickvals = get(gca, 'XTick');
xticklabels = cell(1,numel(xtickvals));
for itick = 1 : numel(xtickvals);
  xticklabels(itick) = cellstr(sprintf('%3.1f',xtickvals(itick)));
end;
set(gca, 'XTickLabel', xticklabels);
xlabel('\beta_Nmax', 'fontsize', 12);
ylabel('Number of shots', 'fontsize', 12);
title('Histograms of \beta_Nmax for disruptive and all shots', 'fontsize', 14);

figure(4);
set(gcf, 'menubar', 'none');
plot(xbins, ybins_disrupt./ybins_all, 's');
ylim([0,1]);
xtickvals = get(gca, 'XTick');
xticklabels = cell(1,numel(xtickvals));
for itick = 1 : numel(xtickvals);
  xticklabels(itick) = cellstr(sprintf('%3.1f',xtickvals(itick)));
end;
set(gca, 'XTickLabel', xticklabels);
ytickvals = get(gca, 'YTick');
yticklabels = cell(1,numel(ytickvals));
for itick = 1 : numel(ytickvals);
  yticklabels(itick) = cellstr(sprintf('%3.1f',ytickvals(itick)));
end;
set(gca, 'YTickLabel', yticklabels);
xlabel('\beta_Nmax', 'fontsize', 12);
ylabel('Disruption fraction', 'fontsize', 12);
title('Fraction of shots that disrupted vs \beta_Nmax', 'fontsize', 14);

index_valid_range = find(betan_disrupt >= 0 & betan_disrupt <= 6);
figure(5);
set(gcf, 'menubar', 'none');
ybins_disrupt2 = hist(betan_disrupt(index_valid_range), xbins);
bar(xbins, ybins_disrupt2, 'FaceColor', 'green');
ylim(ylimits);
xtickvals = get(gca, 'XTick');
xticklabels = cell(1,numel(xtickvals));
for itick = 1 : numel(xtickvals);
  xticklabels(itick) = cellstr(sprintf('%3.1f',xtickvals(itick)));
end;
set(gca, 'XTickLabel', xticklabels);
xlabel('\beta_N at time of disruption', 'fontsize', 12);
ylabel('Number of shots', 'fontsize', 12);
title('Histogram of \beta_N at time of disruption', 'fontsize', 14);

figure(6);
set(gcf, 'menubar', 'none');
bar(xbins, ybins_disrupt, 'r', 'EdgeColor', 'black');
hold on;
bar(xbins+0.05, ybins_disrupt2, 'g', 'EdgeColor', 'black');
hold off;
xlim([0,6]);
ylim(ylimits);
xtickvals = get(gca, 'XTick');
xticklabels = cell(1,numel(xtickvals));
for itick = 1 : numel(xtickvals);
  xticklabels(itick) = cellstr(sprintf('%3.1f',xtickvals(itick)));
end;
set(gca, 'XTickLabel', xticklabels);
xlabel('\beta_N');
ylabel('Number of shots', 'fontsize', 12);
title('Histograms of \beta_Nmax and \beta_N at time of disruption');
