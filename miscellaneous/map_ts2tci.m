function [t, z, n_e, n_e_sig] = map_ts2tci(shot, nlnum)
%   core_mult=core_mult, edge_mult=edge_mult

%if not keyword_set(core_mult) then core_mult=1.0
%if not keyword_set(edge_mult) then edge_mult=1.0

core_mult = 1.0;
edge_mult = 1.0;
%{
map TS density profile onto a given TCI chord (1--10)
INPUTS: shot
;        nlnum = 4 for nl04, etc.
;OUTPUTS: t times
;         z z-mapping onto designated TCI chord
;         n_e and n_e_sig: density and density error bar
%}
t = [1e32];
z = [1e32];
n_e = [1e32];
n_e_sig = [1e32];

[good_indices, efit_times] = efit_check(shot);
if (isempty(good_indices))
  fprintf(1, 'ERROR: No valid EFITs available on shot %i\n', shot);
  return;
end
t1 = min(efit_times);
t2 = max(efit_times);
mdsopen('analysis', shot);
psia = mdsvalue('\efit_aeqdsk:SIBDRY');
psia_t = mdsvalue('dim_of(\efit_aeqdsk:SIBDRY)');
psi0 = mdsvalue('\efit_aeqdsk:SIMAGX');
mdsclose;

[~, status] = mdsopen('electrons', shot);
if (mod(status,2) == 0)
  fprintf(1, 'ERROR: No ELECTRONS tree on shot %i\n', shot);
  return;
end
[nets_core, status] = mdsvalue('.YAG_NEW.RESULTS.PROFILES:NE_RZ');
if (mod(status,2) == 0)
  fprintf(1, 'ERROR: No TS data on shot %i\n', shot);
  mdsclose;
  return;
end
nets_core_err = mdsvalue('.YAG_NEW.RESULTS.PROFILES:NE_ERR');
tts = mdsvalue('dim_of(.YAG_NEW.RESULTS.PROFILES:NE_RZ)');
zts_core = mdsvalue('.YAG_NEW.RESULTS.PROFILES:Z_SORTED');
mts_core = length(zts_core);
zts_edge = mdsvalue('\fiber_z');
mts_edge = length(zts_edge);
[nets_edge, status] = mdsvalue('\ts_ne');
if (mod(status,2) == 0)
  nets_edge = zeros(length(nets_core(:,1)), mts_edge);
  nets_edge_err = nets_edge + 1e20;
else
  nets_edge_err = mdsvalue('\ts_ne_err');
end
mts = mts_core + mts_edge;
rts = mdsvalue('.YAG.RESULTS.PARAM:R') + zeros(1,mts);
[rtci, status] = mdsvalue('.tci.results:rad');
if (mod(status,2) == 0)
  fprintf(1, 'ERROR: No TCI data available on shot %i\n', shot);
  mdsclose;
  return;
end
mdsclose;

nts = length(tts);
zts = zeros(1, mts);
zts(1:mts_core) = zts_core;
zts(mts_core+1:end) = zts_edge;
nets = zeros(nts,mts);
nets_err = zeros(nts,mts);
nets(:, 1:mts_core) = nets_core*core_mult;
nets_err(:, 1:mts_core) = nets_core_err*core_mult;
nets(: ,mts_core+1:end) = nets_edge*edge_mult;
nets_err(:, mts_core+1:end) = nets_edge_err*edge_mult;

ind = find(tts >= t1 & tts <= t2);
if isempty(ind)
  fprintf(1, 'ERROR: Insufficient EFIT times on shot %i', shot);
  return;
end
tts = tts(ind);
nets = nets(ind,:);
nets_err = nets_err(ind,:);

psits=transpose(efit_rz2psi(rts,zts,tts,shot));

mtci = 101;
ztci = -0.4+0.8*[0:mtci-1]/(mtci-1);
rtci = rtci(nlnum)+zeros(1,mtci);

psitci=transpose(efit_rz2psi(rtci,ztci,tts,shot));

psia=interp1(psia_t,psia,tts);
psi0=interp1(psia_t,psi0,tts);
nts = length(tts);
%convert to normalized Psi:
for i=1:nts
  psits(i,:)=(psits(i,:)-psi0(i))/(psia(i)-psi0(i));
  psitci(i,:)=(psitci(i,:)-psi0(i))/(psia(i)-psi0(i));
end  


zmapped=zeros(nts,2*mts)+1e32;
nemapped=zmapped;
nemapped_err=zmapped;

for i=1:nts 
  [psimin,imin]=min(psitci(i,:));
  for j=1:mts 
    if psits(i,j) >= psimin 
      a1=interp1(psitci(i,1:imin),ztci(1:imin),psits(i,j));
      a2=interp1(psitci(i,imin+1:end),ztci(imin+1:end),psits(i,j));
      zmapped(i,[j,j+mts])=[a1,a2];
      nemapped(i,[j,j+mts])=nets(i,j);
      nemapped_err(i,[j,j+mts])=nets_err(i,j);
    end
  end
[~,ind]=sort(zmapped(i,:));
zmapped(i,:)=zmapped(i,ind);
nemapped(i,:)=nemapped(i,ind);
nemapped_err(i,:)=nemapped_err(i,ind);
end

z = zmapped;
n_e = nemapped;
n_e_sig = nemapped_err;
t = tts;

end
  

