if (exist('time_until_disrupt','var')==0 || exist('time','var')==0 || ...
    exist('dipprog_dt','var')==0);
  fprintf(1,['Call the routine ''retrieve_all_params'' before calling ' ...
             'this routine.\n']);
  return;
end;

ipprog = ip - ip_error;
indices_flattop_1 = find(abs(dipprog_dt) <= 6e4);
indices_flattop_2 = find( abs(ipprog) > 1.e5);
indices_flattop = intersect(indices_flattop_1, indices_flattop_2);

fprintf(1,'Use ''indices_flattop'' to select data during flattop.\n');
