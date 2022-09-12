addpath(['/fusion/projects/disruption_warning' ...
    '/matlab_programs']);

% Retrieve DIII-D data and define flattop indices

db = set_database('d3drdb');
%retrieve_all_params_d3d_montes;
retrieve_all_params;
define_indices;

% Subselect shots from the 2015 campaign

blessed_shots=dlmread(['/fusion/projects/disruption_warning/' ...
    'matlab_programs/shotlist_rea_blessed.txt']);
blessed_shots = int32(blessed_shots);
%blessed_shots = int32(dlmread('/home/granetz/Collaborations/DIII-D/disruption_warning_database/shotlist_rea_blessed.txt'));
indices_2015 = find(ismember(shot,blessed_shots));

% Plot li for all flattop disruption times

ii = intersect(indices_flattop,indices_2015);

plot(-time_until_disrupt(ii),li(ii),'.-')
xlabel('Time until disrupt [s]','fontsize', 14)
ylabel('l_i','fontsize', 14)
title('l_i for all flattop disruption times','fontsize', 16)
xlim([-1,0])
ylim([0.8,1.6])