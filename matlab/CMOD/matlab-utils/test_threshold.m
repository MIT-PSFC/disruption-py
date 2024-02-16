% Calculate total number of shots with data in flattop and the total number
% of shots disrupted in the flattop

total_flattop_shots = size(unique(shot(indices_flattop_no_disrupt)),1) + ...
    size(unique(shot(indices_flattop_disrupt_in_flattop)),1);
total_flattop_disruptions = size(unique(shot(...
    indices_flattop_disrupt_in_flattop)),1);
time_before_disrupt = 0.02;

% Specify range of thresholds and calculate false positive rate and
% disruption detection rate at each threshold. 

threshold = [0.8:0.01:2];
false_pos = NaN(size(threshold));
disruptions_detected = NaN(size(threshold));

for i = 1:size(threshold,2);
    disruptions_below_thresh = size(unique(shot(...
        indices_flattop_disrupt_in_flattop(find(li(...
        indices_flattop_disrupt_in_flattop) <= threshold(i) & ...
        time_until_disrupt(indices_flattop_disrupt_in_flattop) ...
        > time_before_disrupt)))),1);
    non_disruptions_below_thresh = size(unique(shot(...
        indices_flattop_no_disrupt(find(li(indices_flattop_no_disrupt) ...
        <= threshold(i))))),1);

    false_pos(i) = non_disruptions_below_thresh/total_flattop_shots;
    disruptions_detected(i) = ...
        disruptions_below_thresh/total_flattop_disruptions;
end

% Plot the false positive and disruption detection rates as a function of
% the threshold chosen and label accordingly

plot(threshold, false_pos, threshold, disruptions_detected)
xlabel('Minimum L_{I} Threshold [H]');
legend('False Positive Rate',['Disruptions detected' ...
    ' 20ms before'],'location','best');