function [p_rad_RT, p_ecrh_RT, p_lh_RT, p_oh_RT, p_icrf_RT, p_nbi_RT,  ...
  rad_input_frac_RT, rad_loss_frac_RT, ip_error_RT,  ...
  q95_RT, beta_p_RT, li_RT, Wmhd_RT] = get_PCS_EAST(shot, timebase);

% This function gets many of the real time signals that are actually used
% in the plasma control system (PCS) on the EAST tokamak.  For development
% of a disruption prediction algorithm that we want to run in real time in
% the PCS, it is better to train on a database of these real time signals
% than on the processed signals in the analysis trees and east trees.
% Note: since 2018 the EFIT-derived signals used in the PCS are calculated
% by P-EFIT, which is different than RT-EFIT.  This information came from
% QP Yuan <qpyuan@ipp.ac.cn>

% Inputs:
%   shot = shot number
%   timebase = times to calculate the powers 
% Outputs:
%   p_rad_RT  = radiated power [W]
%   p_ecrh_RT = electron cyclotron resonance heating power [W]
%   p_lh_RT   = lower hybrid power [W]
%   p_oh_RT   = ohmic power [W]
%   p_icrf_RT = ion cyclotron power [W]
%   p_nbi_RT = neutral beam injection power [W]
%   rad_input_frac = p_RAD/p_input
%                  = p_RAD/(p_OHM + p_LH + p_ICRF + p_ECRH + p_NBI)
%                  = ratio of radiated power to total input power
%   rad_loss_frac  = p_RAD/p_loss
%                  = p_RAD/(p_RAD + p_cond + p_conv)
%                  = p_RAD/(p_input - dWmhd/dt)
%                  = ratio of radiated power to total loss power
%   ip_error_RT    = error between actual and pre-programmed plasma currents
%   q95_RT         = q95 calculated by RT-EFIT for PCS (P-EFIT since 2018)
%   beta_p_RT      = beta_p calculated by RT-EFIT for PCS (P-EFIT since 2018)
%   li_RT          = li calculated by RT-EFIT for PCS (P-EFIT since 2018)
%   Wmhd_RT        = Wmhd calculated by RT-EFIT for PCS (P-EFIT since 2018)
%
% Author: R. Granetz  2018-12-11
% Revision history:
%                          

% We want the output vectors to have the same shape as the input vector,
% "timebase", i.e. if "timebase" is a row vector, the output vectors will
% be row vectors; if "timebase" is a column vector, the output vectors will
% be column vectors.  Note that Matlab calls to mdsvalue return column
% vectors for 1-D signals, so this routine works internally with column
% vectors.  So the first step is to create a column vector copy of the
% input parameter, "timebase".

if (size(timebase,1) > 1);
  timebase_column = timebase;
else;
  timebase_column = transpose(timebase);
end;

p_rad_RT = NaN(size(timebase));
p_ecrh_RT = NaN(size(timebase));
p_lh_RT = NaN(size(timebase));
p_oh_RT = NaN(size(timebase));
p_icrf_RT = NaN(size(timebase));
p_nbi_RT = NaN(size(timebase));
rad_input_frac_RT = NaN(size(timebase));
rad_loss_frac_RT = NaN(size(timebase));
ip_error_RT = NaN(size(timebase));
q95_RT = NaN(size(timebase));
beta_p_RT = NaN(size(timebase));
li_RT = NaN(size(timebase));
Wmhd_RT = NaN(size(timebase));

[shotopened, status] = mdsopen('pcs_east', double(shot));

if (mod(status,2) == 0);
  return;
end;

[timearray, status] = mdsvalue('dim_of(\pcprad)');
if (mod(status,2) == 1);
  signal = mdsvalue('\pcprad');
  p_rad_RT = interp1(timearray, signal, timebase_column, 'linear', 0);
end;

[timearray, status] = mdsvalue('dim_of(\lmeip)');
if (mod(status,2) == 1);
  signal = mdsvalue('\lmeip');
  ip_error_RT = interp1(timearray, signal, timebase_column, 'linear', 0);
end;

[timearray, status] = mdsvalue('dim_of(\pfsbetap)');
if (mod(status,2) == 1);
  signal = mdsvalue('\pfsbetap');
  beta_p_RT = interp1(timearray, signal, timebase_column, 'linear', 0);
end;

[timearray, status] = mdsvalue('dim_of(\pfsli)');
if (mod(status,2) == 1);
  signal = mdsvalue('\pfsli');
  li_RT = interp1(timearray, signal, timebase_column, 'linear', 0);
end;

[timearray, status] = mdsvalue('dim_of(\pfswmhd)');
if (mod(status,2) == 1);
  signal = mdsvalue('\pfswmhd');
  dWmhd_dt = gradient(signal, timearray);
  Wmhd_RT = interp1(timearray, signal, timebase_column, 'linear', 0);
  dWmhd_dt = interp1(timearray, dWmhd_dt, timebase_column, 'linear', 0);
end;

[~, status] = mdsopen('pefitrt_east', shot);
if (mod(status,2) == 1);
  [timearray, status] = mdsvalue('dim_of(\q95)');
  [timearray, indices_unique] = unique(timearray);  % deal with bug
  if (mod(status,2) == 1 && length(timearray) >= 2);
    signal = mdsvalue('\q95');
    signal = signal(indices_unique);                % deal with bug
    q95_RT = interp1(timearray, signal, timebase_column, 'linear', 0);
  end;
  mdsclose;
else;
end;

[timearray, status_I] = mdsvalue('dim_of(\pcnbi1li)');
if (mod(status_I,2) == 1);
  signal = mdsvalue('\pcnbi1li');
  I_nbiL_RT = interp1(timearray, signal, timebase_column, 'linear', 0);
end;
[timearray, status_V] = mdsvalue('dim_of(\pcnbi1lv)');
if (mod(status_V,2) == 1);
  signal = mdsvalue('\pcnbi1lv');
  V_nbiL_RT = interp1(timearray, signal, timebase_column, 'linear', 0);
end;
if (mod(status_I,2) == 1 && mod(status_V,2) == 1);
  p_nbiL_RT = I_nbiL_RT .* V_nbiL_RT;
end;

[timearray, status_I] = mdsvalue('dim_of(\pcnbi1ri)');
if (mod(status_I,2) == 1);
  signal = mdsvalue('\pcnbi1ri');
  I_nbiR_RT = interp1(timearray, signal, timebase_column, 'linear', 0);
end;
[timearray, status_V] = mdsvalue('dim_of(\pcnbi1rv)');
if (mod(status_V,2) == 1);
  signal = mdsvalue('\pcnbi1rv');
  V_nbiR_RT = interp1(timearray, signal, timebase_column, 'linear', 0);
end;
if (mod(status_I,2) == 1 && mod(status_V,2) == 1);
  p_nbiR_RT = I_nbiR_RT .* V_nbiR_RT;
  p_nbi_RT = p_nbiL_RT + p_nbiR_RT;
end;

[timearray, status_I] = mdsvalue('dim_of(\pcplhi)');
if (mod(status_I,2) == 1);
  signal = mdsvalue('\pcplhi');
  p_lh_46_inj_RT = interp1(timearray, signal, timebase_column, 'linear', 0);
end;
[timearray, status_R] = mdsvalue('dim_of(\pcplhr)');
if (mod(status_R,2) == 1);
  signal = mdsvalue('\pcplhr');
  p_lh_46_ref_RT = interp1(timearray, signal, timebase_column, 'linear', 0);
end;
if (mod(status_I,2) == 1 && mod(status_R,2) == 1);
  p_lh_46_RT = p_lh_46_inj_RT - p_lh_46_ref_RT;
end;

[timearray, status_I] = mdsvalue('dim_of(\pcplhi2)');
if (mod(status_I,2) == 1);
  signal = mdsvalue('\pcplhi2');
  p_lh_245_inj_RT = interp1(timearray, signal, timebase_column, 'linear', 0);
end;
[timearray, status_V] = mdsvalue('dim_of(\pcplhr2)');
if (mod(status_V,2) == 1);
  signal = mdsvalue('\pcplhr2');
  p_lh_245_ref_RT = interp1(timearray, signal, timebase_column, 'linear', 0);
end;
if (mod(status_I,2) == 1 && mod(status_R,2) == 1);
  p_lh_245_RT = p_lh_245_inj_RT - p_lh_245_ref_RT;
  p_lh_RT = p_lh_46_RT + p_lh_245_RT;
end;

%{
 Q.P. Yuan:
There is no signals from ICRF or ECRH connected to PCS. And I checked the
signals for NBI, there is no effective value, just noise. We will make
the NBI signals available in this EAST campaign. The PLHI2 and PLHR2 for
2.45G  LHW system have signals but are not calibrated yet. 
%}

P_ohm = get_P_ohm_EAST(shot, timebase);  % Note, I'm not bothering to
                                         % to calculate a real time version
mdsclose;

%{
% Total input power [W]
p_input_RT = P_ohm + p_lh_RT + p_icrf_RT + p_ecrh_RT + p_nbi_RT;
rad_input_frac_RT = p_rad_RT ./ p_input_RT;
rad_loss_frac_RT  = p_rad_RT ./ (p_input_RT - dWmhd_dt);
%}

size_timebase = size(timebase);
if (size(p_rad_RT) ~= size_timebase);
  p_rad_RT = transpose(p_rad_RT);
end;
if (size(p_ecrh_RT) ~= size_timebase);
  p_ecrh_RT = transpose(p_ecrh_RT);
end;
if (size(p_lh_RT) ~= size_timebase);
  p_lh_RT = transpose(p_lh_RT);
end;
if (size(p_oh_RT) ~= size_timebase);
  p_oh_RT = transpose(p_oh_RT);
end;
if (size(p_icrf_RT) ~= size_timebase);
  p_icrf_RT = transpose(p_icrf_RT);
end;
if (size(p_nbi_RT) ~= size_timebase);
  p_nbi_RT = transpose(p_nbi_RT);
end;
if (size(rad_input_frac_RT) ~= size_timebase);
  rad_input_frac_RT = transpose(rad_input_frac_RT);
end;
if (size(rad_loss_frac_RT) ~= size_timebase);
  rad_loss_frac_RT = transpose(rad_loss_frac_RT);
end;
if (size(ip_error_RT) ~= size_timebase);
  ip_error_RT = transpose(ip_error_RT);
end;
if (size(q95_RT) ~= size_timebase);
  q95_RT = transpose(q95_RT);
end;
if (size(beta_p_RT) ~= size_timebase);
  beta_p_RT = transpose(beta_p_RT);
end;
if (size(li_RT) ~= size_timebase);
  li_RT = transpose(li_RT);
end;
if (size(Wmhd_RT) ~= size_timebase);
  Wmhd_RT = transpose(Wmhd_RT);
end;

end
