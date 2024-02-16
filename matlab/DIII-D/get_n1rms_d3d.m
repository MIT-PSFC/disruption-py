function [n1rms, n1rms_normalized] = get_n1rms_d3d(shot, timebase)
if (size(timebase,2) > 1)
  timebase = transpose(timebase);
end

n1rms = NaN(size(timebase));
n1rms_normalized = NaN(size(timebase));

mdsconnect('atlas.gat.com');
[~, status] = mdsopen('d3d', shot);
if (mod(status,2) == 0)
  return
end

[time, status] = mdsvalue('dim_of(\n1rms)');
time = time/1.e3; % convert ms to s

if (mod(status,2)==0 || length(time) <= 4)
  return
end
n1rms = mdsvalue('\n1rms');
n1rms = n1rms * 1e-4;  % convert gauss to tesla
n1rms = interp1(time, n1rms, timebase, 'linear');

[Btime, status] = mdsvalue(['dim_of(ptdata("bt", ' num2str(shot) '))']);
if (mod(status,2) == 1)
  Btime = Btime/1.e3; % convert ms to s
  Btor = mdsvalue(['ptdata("bt", ' num2str(shot) ')']); %[T]
  mdsclose;
  Btor = interp1(Btime, Btor, timebase, 'linear');
  n1rms_normalized = n1rms ./ abs(Btor);
else
  mdsclose;  
  return;
end

end
