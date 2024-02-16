function [indices_thinned] = thin_points(time_until_disrupt, start_time, ...
  remainder_fraction);

indices_to_be_thinned = find(time_until_disrupt <= start_time & ...
  time_until_disrupt >= 0);
indices_not_to_be_thinned = find(time_until_disrupt > start_time);

initial_number_of_indices_to_be_thinned = length(indices_to_be_thinned);
desired_number_of_thinned_indices = ...
  round(initial_number_of_indices_to_be_thinned * remainder_fraction);

% if row vector input, then row vector output, and vice versa

if size(time_until_disrupt, 2) > 1;
  indices_thinned = [indices_not_to_be_thinned, ...
    randsample(indices_to_be_thinned, desired_number_of_thinned_indices)];
else;
  indices_thinned = [indices_not_to_be_thinned; ...
    randsample(indices_to_be_thinned, desired_number_of_thinned_indices)];
end;
