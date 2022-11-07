function [nlts, nlts_t] = integrate_ts2tci(shot, nlnum);
%  graph=graph, core_mult=core_mult, edge_mult=edge_mult

%if not keyword_set(core_mult) then core_mult=1.0
%if not keyword_set(edge_mult) then edge_mult=1.0

core_mult = 1.0;
edge_mult = 1.0;

nlts = 1e32;
nlts_t = 1e32;

[t, z, n_e, n_e_sig] = map_ts2tci(shot, nlnum);
%   core_mult=core_mult, edge_mult=edge_mult

if (z(1,1) == 1e32);
  fprintf(1, 'ERROR: No valid TS mappings to TCI chord for shot %i\n', shot);
  return;
end;

nts = length(t);
nlts_t = t;
nlts = NaN(size(t));

for i = 1:nts;
  ind = find(abs(z(i,:)) < 0.5 & n_e(i,:) > 0 & n_e(i,:) < 1e21 & ...
            n_e(i,:) ./ n_e_sig(i,:) > 2);
  if (length(ind) < 3);
    nlts(i)=0.;
  else;
    x = z(i,ind);
    y = n_e(i,ind);
    ind_uniq = unique(x);
    [x, ind_uniq, ~] = unique(x);
    y = y(ind_uniq);
    nlts(i) = trapz(x,y);
  end;
% plot,z(i,*),n_e(i,*),xrange=[-.5,.5]
% oploterr,z(i,*),n_e(i,*),n_e_sig(i,*)
% oplot,z(i,ind),n_e(i,ind),psym=4,col=2
% stop
end;

%{
node = ['.TCI.RESULTS:NL_' num2str(nlnum, '%02i')];
plot(nlts_t, nlts, '^');
xlabel('Time (s)');
ylabel('\int n_e dL (10^{20} m^{-3})');
title = (node);
[tci, tci_t] = mdsget('electrons',shot,node);
hold on;
plot(tci_t, tci, 'r');
%}

end
