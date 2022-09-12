% Retrieve C-Mod data and define flattop indices

db = set_database('logbook');
retrieve_all_params;
define_indices;

% Subselect shots from the 2015 campaign

indices_2015 = find(shot>1150000000 & shot<1160000000);

% Plot li for all flattop disruption times

ii = intersect(indices_flattop,indices_2015);

plot(-time_until_disrupt(ii),li(ii),'.')
xlabel('Time until disrupt [s]','fontsize', 14)
ylabel('l_i','fontsize', 14)
title('l_i for all flattop disruption times','fontsize', 16)
xlim([-1,0])
ylim([0.8,1.6])