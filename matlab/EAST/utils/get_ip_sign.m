if (exist('shot','var')==0 || exist('ip','var')==0);
  fprintf(1,['Call the routine ''retrieve_all_params'' before calling ' ...
             'this routine.\n']);
  return;
end;

shotlist = unique(shot);
ipsign = NaN(size(ip));
for ishot = 1:length(shotlist);
  qq = find(shot == shotlist(ishot));
  plasma_current = ip(qq);
  plasma_current = plasma_current(find(~isnan(plasma_current)));
  ip_integrated = sum(plasma_current);
  ipsign(qq) = ones(size(qq)) * sign(ip_integrated);
end;

fprintf(1,'The variable ''ipsign'' contains the algebraic sign of Ip.\n');
