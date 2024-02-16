function [signal, dim_of_signal] = mdsget(tree, shot, node);

% Translated from IDL to Matlab by R. Granetz -- 2019/01/17

% get the value of a node in a given shot with one line.

signal = [];
dim_of_signal =[];

[~, status] = mdsopen(tree, shot);
if (mod(status,2) == 1);
  [result, status] = mdsvalue(node);
  if (mod(status,2) == 0);
    fprintf(1, 'ERROR: Unable to read data from node ''%s''\n', node);
  else;
    signal = result;
    [dim_of_signal, status] = mdsvalue(['dim_of(' node ')']);
  end;
  mdsclose;
else;
  fprintf(1, 'ERROR: Unable to open tree ''%s'' on shot %i\n', tree, shot);
end;

end
