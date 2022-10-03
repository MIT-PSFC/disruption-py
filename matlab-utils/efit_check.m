function [good_indices, efit_times] = efit_check(shot);

[~, status] = mdsopen('analysis', shot);
if (mod(status,2) == 0);
  good_indices = [];
  efit_times = [];
  return;
end;

test = ['_lf=\analysis::efit_aeqdsk:lflag,' ...
        '_l0=((sum(_lf,1) - _lf[*,20] - _lf[*,1])==0),' ...
        '_n=\analysis::efit_fitout:nitera,(_l0 and (_n>4))'];

[dummy, status] = mdsvalue(test);
[good_indices, ntimes] = find(dummy);

if (mod(status,2) == 1 & ntimes > 0);
  tt = mdsvalue('dim_of(_lf)');
  efit_times = tt(good_indices);
else;
  good_indices = [];
  efit_times = [];
end;
mdsclose;

end
