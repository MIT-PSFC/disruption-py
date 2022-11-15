function [ip_new, q0_new, rad_frac_new, ml_new, li_new] = get_vars(ip, q0, rad_frac, ml, li, time, timebase);
  ip_new = interp1(time, ip, timebase, 'linear');
  q0_new = interp1(time, q0, timebase, 'linear');
  rad_frac_new = interp1(time, rad_frac, timebase, 'linear');
  ml_new = interp1(time, ml, timebase, 'linear');
  li_new = interp1(time, li, timebase, 'linear');

end

