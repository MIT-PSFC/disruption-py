function psi = efit_rz2psi(r, z, t, shot, tree);

psi = NaN(length(r), length(t));
z=single(z);
if (nargin == 5);
  [~, status] = mdsopen(tree, shot);
  if (mod(status,2) == 0);
    fprintf(1, 'ERROR: no %s tree for shot %i\n', tree, shot);
    return;
  end;
elseif (nargin == 4);
  tree = 'analysis';
  [~, status] = mdsopen(tree, shot);
  if (mod(status,2) == 0);
    fprintf(1, 'ERROR: no %s tree for shot %i\n', tree, shot);
    return;
  end;
elseif (nargin == 3);

% Assume the proper EFIT tree is already open

elseif (nargin < 3);
  fprintf(1, ['ERROR: at least the first 3 arguments must be\n' ...
              '       specified when calling EFIT_RZ2PSI\n']);
  psi = '';
  return;
end;
[~, status] = mdsvalue('\efit_geqdsk:psirz');
%[psirz, status] = 2*pi*mdsvalue('\efit_geqdsk:psirz');
if (mod(status,2) == 0);
  fprintf(1, 'ERROR: no EFIT psirz data for specified shot %i\n');
  mdsclose;
  return;
else;
  psirz = 2*pi*mdsvalue('\efit_geqdsk:psirz');
  rgrid = mdsvalue('dim_of(\efit_geqdsk:psirz,0)');
  zgrid = mdsvalue('dim_of(\efit_geqdsk:psirz,1)');
  times = mdsvalue('dim_of(\efit_geqdsk:psirz,2)');
  mdsclose;
end;

[Rgrid, Zgrid] = meshgrid(rgrid, zgrid);

for itime = 1:length(t);

% Select EFIT times closest to the requested time

[~, indx] = min(abs(times - t(itime)));
Psirz = transpose(psirz(:,:,indx));

psi(:, itime) = interp2(Rgrid, Zgrid, Psirz, r, z, 'cubic');

end;

end
