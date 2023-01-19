function fit_TS_polynomial(shot, times);

mdsconnect('atlas.gat.com');

[shotopened, status]=mdsopen('electrons', shot);
if (mod(status,2)==0);
% fprintf(1,'  Unable to open ELECTRONS tree for shot%7i\n', shot);
  return;
end;

% Read in Thomson core temperature data, which is a 2-D array, with the
% dependent dimensions being time and z (vertical coordinate)

[TS_time, status] = mdsvalue('dim_of(\top.ts.blessed.core:temp,0)');

if (mod(status,2) == 1);
  TS_data = mdsvalue('\top.ts.blessed.core:temp');
  TS_z = mdsvalue('dim_of(\top.ts.blessed.core:temp,1)');
  TS_time = TS_time/1.e3; % convert ms to s
% Get rid of the last channel (#41), which is not real
  TS_data = TS_data(:, 1:end-1);
  TS_z = TS_z(1:end-1);
  mdsclose;
else;
  mdsclose;
  return;
end;

xarray = [0:.05:.9];
itimes = find(TS_time > 0);
yarray = NaN(length(itimes), length(xarray));

% plot([0,1],[0,1],'w'); xlim([0,.9]); ylim([0,3000]); hold on;

for i = 1:length(itimes);
  y = TS_data(itimes(i),:);
  ok_indices = find(y ~= 0);
  y = y(ok_indices);
  x = TS_z(ok_indices);
  if (length(ok_indices) > 2);
    p = polyfit(x, transpose(y), 2);
    yarray(i,:) = polyval(p, xarray);
%   plot(TS_z, TS_data(itimes(i),:), 'sb');
%   plot(xarray, yarray(i,:), 'r');
%   fprintf(1, '%10.3f\n', TS_time(itimes(i)));
%   pause(0.1);
%   plot(TS_z, TS_data(itimes(i),:), 'sw');
%   plot(xarray, yarray(i,:), 'w');
  end;
end;

end
