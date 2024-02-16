if (exist('time_until_disrupt','var')==0 || exist('time','var')==0 || ...
    exist('dipprog_dt','var')==0);
  fprintf(1,['Call the routine ''retrieve_all_params'' before calling ' ...
             'this routine.\n']);
  return;
end;

indices_flattop = find(abs(dipprog_dt) <= 1.e3); % <= 1 kA/s limit

fprintf(1,'Use ''indices_flattop'' to select data during flattop.\n');
