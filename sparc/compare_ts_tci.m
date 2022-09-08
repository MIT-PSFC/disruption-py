function [nlts1, nlts2, nltci1, nltci2,time1,time2]=compare_ts_tci(shot, nlnum)

%   core_mult=core_mult, edge_mult=edge_mult

% if not keyword_set(core_mult) then core_mult=1.0
% if not keyword_set(edge_mult) then edge_mult=1.0

core_mult = 1.0;
edge_mult = 1.0;

nlts1 = [1e32];
nlts2 = [1e32];
nltci1 = [1e32];
nltci2 = [1e32];
tstime1 = [1e32];
tstime2 = [1e32];
shotstr=strtrim(int2str(shot));
tcinode = ['\top.electrons.TCI.RESULTS:NL_' num2str(nlnum, '%02i')];

%[~, status] = mdsopen('electrons', shot);
%if (mod(status,2) == 0);
%  fprintf(1, 'ERROR: unable to open ELECTRONS tree for shot %i\n', shot);
%  return
%else;
  [tstime, status] = mdsvalue('dim_of(\top.electrons.YAG_NEW.RESULTS.PROFILES:NE_RZ)');
  if (mod(status,2) == 0);
    fprintf(1, 'ERROR: unable to read Thomson data for shot %i\n', shot);
%   mdsclose;
    return;
  end;
  [tci, status] = mdsvalue(tcinode);
  if (mod(status,2) == 0);
    fprintf(1, shotstr+': No TCI.');
%   mdsclose;
    return;
  end;
  tci_t = mdsvalue(['dim_of(' tcinode ')']);
end;
%mdsclose;
  
[nlts, nlts_t] = integrate_ts2tci(shot, nlnum);
%  core_mult=core_mult, edge_mult=edge_mult
t0 = min(nlts_t);
t1 = max(nlts_t);

[nyag1, nyag2, indices1, indices2] = parse_yags(shot);

if (nyag1 > 0);
    indices1=indices1+1;
  tstime1 = tstime(indices1);
  ind = find(tstime1 >= t0 & tstime1 <= t1);
  if (length(ind) > 0);
    nltci1 = interp1(tci_t, tci, tstime1(ind));
    nlts1 = interp1(nlts_t, nlts, tstime1(ind));
    time1 = tstime1(ind);
  end;
else;
  time1 = -1;
end;

if (nyag2 > 0);
    indices2=indices2+1;
  tstime2 = tstime(indices2);
  ind = find(tstime2 >= t0 & tstime2 <= t1);
  if (length(ind) > 0);
    nltci2 = interp1(tci_t, tci, tstime2(ind));
    nlts2 = interp1(nlts_t, nlts, tstime2(ind));
    time2 = tstime2(ind);
  end;
else;
  time2=-1;
end;

end
