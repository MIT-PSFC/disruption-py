function [upper_gap, lower_gap] = get_PEFIT_gaps(shot)

mdsconnect('mds.ipp.ac.cn');

[~, status] = mdsopen('pefit_east', double(shot));
if (mod(status,2) == 0);
  [~, status] = mdsopen('efit_east', double(shot));
  if (mod(status,2) == 0);
    error('EFIT tree error');
  end;
end;

[efittime, status] = mdsvalue('dim_of(\top.results.geqdsk:bdry)');
[efittime, ~, ~] = unique(efittime);
if (mod(status,2)==0 || length(efittime) <= 1);
  error('EFIT timebase error');
end;
upper_gap  = NaN(size(efittime));
lower_gap  = NaN(size(efittime));

[data, status] = mdsvalue('\top.results.geqdsk:bdry');
if (mod(status,2)==0);
  return
end;
xcoords = reshape(data(1,:,:), [], length(efittime));
ycoords = reshape(data(2,:,:), [], length(efittime));

[xfirstwall, status] = mdsvalue('\top.results.geqdsk:xlim');
if (mod(status,2)==0);
  return
end;
yfirstwall = mdsvalue('\top.results.geqdsk:ylim');
mdsclose();

seed = ones(size(xcoords,1),1);
xfirstwall = reshape(xfirstwall,[],1);
yfirstwall = reshape(yfirstwall,[],1);
xfirstwall_mat = repmat(xfirstwall,1,length(efittime));
yfirstwall_mat = repmat(yfirstwall,1,length(efittime));

index_upperwall = find(yfirstwall>0.6);
index_lowerwall = find(yfirstwall<-0.6);
xupperwall = xfirstwall_mat(index_upperwall,:);
xupperwall_mat = reshape(kron(xupperwall,seed),length(seed),[],length(efittime));
xlowerwall = xfirstwall_mat(index_lowerwall,:);
xlowerwall_mat = reshape(kron(xlowerwall,seed),length(seed),[],length(efittime));
yupperwall = yfirstwall_mat(index_upperwall,:);
yupperwall_mat = reshape(kron(yupperwall,seed),length(seed),[],length(efittime));
ylowerwall = yfirstwall_mat(index_lowerwall,:);
ylowerwall_mat = reshape(kron(ylowerwall,seed),length(seed),[],length(efittime));
%index_upperplasma = find(ycoords(:,1)>0);
%index_lowerplasma = find(ycoords(:,1)<0);
%xupperplasma = xcoords(index_upperplasma,:);
%xlowerplasma = xcoords(index_lowerplasma,:);
%yupperplasma = ycoords(index_upperplasma);
%ylowerplasma = ycoords(index_lowerplasma);
xupperplasma_mat = reshape(repmat(xcoords,size(xupperwall,1),1),[],size(xupperwall,1),length(efittime));
yupperplasma_mat = reshape(repmat(ycoords,size(yupperwall,1),1),[],size(yupperwall,1),length(efittime));
xlowerplasma_mat = reshape(repmat(xcoords,size(xlowerwall,1),1),[],size(xlowerwall,1),length(efittime));
ylowerplasma_mat = reshape(repmat(ycoords,size(ylowerwall,1),1),[],size(ylowerwall,1),length(efittime));

uppergap_mat = (((xupperplasma_mat-xupperwall_mat).^2) + ((yupperplasma_mat-yupperwall_mat).^2)).^0.5;
lowergap_mat = (((xlowerplasma_mat-xlowerwall_mat).^2) + ((ylowerplasma_mat-ylowerwall_mat).^2)).^0.5;
upper_gap = reshape(min(min(uppergap_mat)),[],1);
lower_gap = reshape(min(min(lowergap_mat)),[],1);

end
