retrieve_all_params;
define_indices;
clearvars -except qstar q95 li kappa indices_* time*;
nd = indices_flattop_no_disrupt;
dd_near = find(time_until_disrupt > 0.005 & time_until_disrupt < 0.35); 
dd_far = find(time_until_disrupt > 0.35);  
tot = length(indices_flattop);                      
tot_nd = length(indices_no_disrupt);
tot_dd = length(indices_disrupt);

bins = [0.:0.05:2.5];
[N_near,edges_near] = histcounts(kappa(dd_near),bins);
[N_far,edges_far] = histcounts(kappa(dd_far),bins);
[N_nd,edges_nd] = histcounts(kappa(nd),bins);

figure();
histogram('BinCounts',N_nd,'BinEdges',edges_nd,'normalization','probability','DisplayStyle', 'stairs', 'EdgeColor', 'g');
xlabel('elongation (\kappa)', 'fontsize', 14);
ylabel('fraction of flattop time slices', 'fontsize', 14);
title('DIII-D', 'fontsize','14');
hold on;
histogram('BinCounts',N_far,'BinEdges',edges_far,'normalization','probability','DisplayStyle', 'stairs', 'EdgeColor', 'b');
hold on;                                                                 
histogram('BinCounts',N_near,'BinEdges',edges_near,'normalization','probability','DisplayStyle', 'stairs', 'EdgeColor', 'r');
legend('non-disruptive', 'time until disrupt > 0.3 s', ...
  '0.005 < time until disrupt <  0.3 s', 'location', 'northwest');
hold off;

%% normalize wrt to all time slices? or just consider flattop time slices?
bins = [0.:0.1:3];
[N_near,edges_near] = histcounts(li(dd_near),bins);
[N_far,edges_far] = histcounts(li(dd_far),bins);
[N_nd,edges_nd] = histcounts(li(nd),bins);

figure();
histogram('BinCounts',N_nd/tot_nd,'BinEdges',edges_nd,'DisplayStyle', 'stairs', 'EdgeColor', 'g');
xlabel('l_i', 'fontsize', 14);
ylabel('fraction of flattop time slices', 'fontsize', 14);
title('DIII-D', 'fontsize','14');
hold on;
histogram('BinCounts',N_far/tot_dd,'BinEdges',edges_far,'DisplayStyle', 'stairs', 'EdgeColor', 'b');
hold on;
histogram('BinCounts',N_near/tot_dd,'BinEdges',edges_near,'DisplayStyle', 'stairs', 'EdgeColor', 'r');
legend('non-disruptive', 'time until disrupt > 0.3 s', ...
  '0.005 < time until disrupt <  0.3 s', 'location', 'northwest');
hold off;

bins = [0.9:0.02:2.03];
figure();
histogram(kappa(nd),bins,'normalization','probability','DisplayStyle', 'stairs', 'EdgeColor', 'g');
xlabel('elongation (\kappa)', 'fontsize', 14);
ylabel('fraction of flattop time slices', 'fontsize', 14);
title('DIII-D', 'fontsize','14');
hold on;
histogram(kappa(dd_far),bins,'normalization','probability','DisplayStyle', 'stairs', 'EdgeColor', 'b');
hold on;                                                          
histogram(kappa(dd_near),bins,'normalization','probability','DisplayStyle', 'stairs', 'EdgeColor', 'r');
legend('non-disruptive', 'time until disrupt > 0.35 s', ...
  '0.005 < time until disrupt <  0.35 s', 'location', 'northwest');
hold off;

bins = [0.4:0.04:2.5];
figure();
histogram(li(nd),bins,'normalization','probability','DisplayStyle', 'stairs', 'EdgeColor', 'g');
xlabel('l_i', 'fontsize', 14);
ylabel('fraction of flattop time slices', 'fontsize', 14);
title('DIII-D', 'fontsize','14');
hold on;
histogram(li(dd_far),bins,'normalization','probability','DisplayStyle', 'stairs', 'EdgeColor', 'b');
hold on;
histogram(li(dd_near),bins,'normalization','probability','DisplayStyle', 'stairs', 'EdgeColor', 'r');
legend('non-disruptive', 'time until disrupt > 0.35 s', ...
  '0.005 < time until disrupt <  0.35 s', 'location', 'northwest');
xlim([0.4,2.5])
hold off;

