function [nyag1, nyag2, indices1, indices2] = parse_yags(shot);

% Translated from IDL to Matlab by R. Granetz -- 2018/11/29

% produces indices used to separate edgets signals from two YAGs, based
% on the number of YAG points in each laser ( nyag1, nyag2 ) and the
% total number of pulses in the edge TS digitizer

mdsopen('electrons',shot);
nyag1 = mdsvalue('\knobs:pulses_q');
[nyag2, status] = mdsvalue('\knobs:pulses_q_2');
if (mod(status,2)==0); nyag2=0; end;
dark = mdsvalue('\n_dark_prior');
ntotal = mdsvalue('\n_total');
nt = ntotal-dark;
mdsclose;

if (nyag1 == 0);
  if (nyag2 == 0);
    indices1 = -1;
    indices2 = -1;
  else;
    indices1 = -1;
    indices2 = [0:nyag2-1];
  end;
else;
  if (nyag2 == 0);
    indices1 = [0:nyag1-1];
    indices2 = -1;
  else;
    if (nyag1 == nyag2);
      indices1 = 2*[0:nyag1-1];
      indices2 = indices1+1;
    else;
      if (nyag1 < nyag2);
        indices1 = 2*[0:nyag1-1];
        indices2 = [2*[0:nyag1-1] + 1, 2*nyag1 + [0:(nyag2-nyag1)-1]];
      else;
        indices2 = 2*[0:nyag2-1]+1;
        indices1 = [2*[0:nyag2-1], 2*nyag2+[0:(nyag1-nyag2)-1]];
      end;
    end;
  end;
end;

ind=find(indices1 < nt);
if (nyag1 > 0) && (length(ind) > 0);
  indices1 = indices1(ind);
else;
  indices1 = -1;
end;
ind = find(indices2 < nt);
if (nyag2 > 0) && (length(ind) > 0);
  indices2 = indices2(ind);
else;
  indices2 = -1;
end;

end
