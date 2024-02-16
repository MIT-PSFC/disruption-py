function diperrdt = iperr_deriv(shots, times, iperr);

% This routine will output an array, 'diperrdt', containing the time
% derivatives of the iperr data.  The output array will have the same
% ordering as the 'shots', 'times', and 'iperr' arrays.

diperrdt = zeros(size(iperr));
shotlist = unique(shots);
for ishot = 1:length(shotlist);
  indices = find(shots == shotlist(ishot));
  diperrdt(indices) = gradient(iperr(indices), times(indices));
end;
infvals = find(isinf(diperrdt));
diperrdt(infvals) = NaN;

end
