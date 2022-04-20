function H98_y2 = get_H98_EAST(shot,timebase)

if (size(timebase,2)>1)
    timebase=transpose(timebase);
end
H98_y2=NaN(size(timebase));

mdsconnect('mds.ipp.ac.cn');
[~, status] = mdsopen('energy_east', double(shot));
if (mod(status,2) == 0)
  return
end

[time, status] = mdsvalue('dim_of(\h98_mhd)');
if (mod(status,2)==0 || length(time) <= 1)
    return
end

[H98_y2, status] = mdsvalue('\h98_mhd');
if (mod(status,2)==0 || length(H98_y2) <= 1)
    return
end
H98_y2=interp1(time,H98_y2,timebase,'linear');

mdsclose();

end
