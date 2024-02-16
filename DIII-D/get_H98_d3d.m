function H98_y2 = get_H98_d3d(shot, timebase)
if (size(timebase,2) > 1)
  timebase = transpose(timebase);
end

H98_y2 = NaN(size(timebase));


%mdsconnect('atlas.gat.com');
[~, status] = mdsopen('transport', shot);
if (mod(status,2) == 0)
  return
end

[time, status] = mdsvalue('dim_of(\H_THH98Y2)');
if (mod(status,2)==0 || length(time) <= 2)
  return
end

time = time/1.e3; % convert ms to s

[H98_y2, status] = mdsvalue('\H_THH98Y2');
if (mod(status,2)==0 || length(H98_y2) <= 2)
  H98_y2 = NaN(size(timebase));
  return
end

H98_y2 = interp1(time, H98_y2, timebase, 'linear');

%mdsclose;

end
