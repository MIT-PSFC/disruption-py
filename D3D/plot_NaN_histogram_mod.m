db = set_database('d3drdb');
retrieve_all_params;
define_indices;

shotlist=dlmread('../new_flattop_disruptions_cleaned.txt');
shotlist=int32(shotlist);
indx_clean_shots = find(ismember(shot, shotlist));
indices_flattop_disrupt_in_flattop = ...
  intersect(indices_flattop_disrupt_in_flattop, indx_clean_shots);

indx_good_ip                 = find(~isnan(ip));
indx_good_ip_error           = find(~isnan(ip_error));
indx_good_n_e                = find(~isnan(n_e));
indx_good_Greenwald_fraction = find(~isnan(Greenwald_fraction));
indx_good_beta_p             = find(~isnan(beta_p));
indx_good_n_equal_1_mode     = find(~isnan(n_equal_1_mode));
indx_good_q95                = find(~isnan(q95));
indx_good_li                 = find(~isnan(li));
indx_good_dWmhd_dt           = find(~isnan(dWmhd_dt));
indx_good_p_rad              = find(~isnan(p_rad));
indx_good_v_loop             = find(~isnan(v_loop));

indx_good = indices_flattop_disrupt_in_flattop;

indx_good = intersect( indx_good, indx_good_ip);
indx_good = intersect( indx_good, indx_good_ip_error);
indx_good = intersect( indx_good, indx_good_n_e);
indx_good = intersect( indx_good, indx_good_Greenwald_fraction);
indx_good = intersect( indx_good, indx_good_beta_p);
indx_good = intersect( indx_good, indx_good_n_equal_1_mode);
indx_good = intersect( indx_good, indx_good_q95);
indx_good = intersect( indx_good, indx_good_li);
indx_good = intersect( indx_good, indx_good_dWmhd_dt);
indx_good = intersect( indx_good, indx_good_p_rad);
indx_good = intersect( indx_good, indx_good_v_loop);

bins = [-0.250 : 0.002 : 0.002] - 0.001;

histogram( -time_until_disrupt(indx_good), bins); xlim([-0.1, 0.0]);
set(gca, 'fontsize', 12);
xlabel('Time before disrupt [s]', 'fontsize', 14);
title({'Histograms of # of records containing ', ...
       'no NaN''s (blue) and at least one NaN (brown)', ...
       'for 11 selected parameters', ...
       '(clean flattop disruptions only)'}, ...
       'fontsize',15);

ii = indices_flattop_disrupt_in_flattop;

indx_bad_ip                  = intersect(ii, find(isnan(ip)));
indx_bad_ip_error            = intersect(ii, find(isnan(ip_error)));
indx_bad_n_e                 = intersect(ii, find(isnan(n_e)));
indx_bad_Greenwald_fraction  = intersect(ii, find(isnan(Greenwald_fraction)));
indx_bad_beta_p              = intersect(ii, find(isnan(beta_p)));  
indx_bad_n_equal_1_mode      = intersect(ii, find(isnan(n_equal_1_mode)));
indx_bad_q95                 = intersect(ii, find(isnan(q95)));
indx_bad_li                  = intersect(ii, find(isnan(li)));
indx_bad_dWmhd_dt            = intersect(ii, find(isnan(dWmhd_dt)));
indx_bad_p_rad               = intersect(ii, find(isnan(p_rad)));
indx_bad_v_loop              = intersect(ii, find(isnan(v_loop)));

indx_bad = [];

indx_bad = union( indx_bad , indx_bad_ip);
indx_bad = union( indx_bad , indx_bad_ip_error);
indx_bad = union( indx_bad , indx_bad_n_e);
indx_bad = union( indx_bad , indx_bad_Greenwald_fraction);
indx_bad = union( indx_bad , indx_bad_beta_p);
indx_bad = union( indx_bad , indx_bad_n_equal_1_mode);
indx_bad = union( indx_bad , indx_bad_q95);
indx_bad = union( indx_bad , indx_bad_li);
indx_bad = union( indx_bad , indx_bad_dWmhd_dt);
indx_bad = union( indx_bad , indx_bad_p_rad);
indx_bad = union( indx_bad , indx_bad_v_loop);

hold on;
histogram( -time_until_disrupt(indx_bad ), bins);
hold off;
