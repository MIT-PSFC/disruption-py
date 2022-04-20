function [p_RAD, p_ECRH, p_LH, p_OHM, p_ICRF, p_NBI, rad_input_frac, ...
  rad_loss_frac, p_input] = get_power_EAST(shot, timebase);

% This function gets the input heating powers -- ohmic (p_OHM),electron 
% cyclotron resonance heating (p_ECRH), neutral beam injection system (p_NBI)
% ion cyclotron (p_ICRF), and lower hybrid (p_LH) -- as well as the radiated
% output power (p_RAD).  If any of the auxiliary heating powers are not
% available (there was no ICRF or LH), then this function returns an array
% of zeros for them.  The ohmic heating power is obtained by calling a
% separate function, get_P_ohm.m.  If either the ohmic heating power or the
% radiated power is unavailable, then arrays of NaN (Not-a-Number) are
% returned for them, and for the radiated power fraction.
% 
% Inputs:
%   shot = shot number
%   timebase = times to calculate the powers 
% Outputs:
%   p_RAD  = radiated power [W]
%   p_ECRH = electron cyclotron resonance heating power [W]
%   p_LH   = lower hybrid power [W]
%   p_OHM  = ohmic power [W]
%   p_ICRF = ion cyclotron power [W]
%   p_NBI = neutral beam injection power [W]
%   p_input = total input power = p_OHM + p_LH + p_ICRF + p_ECRH + p_NBI [W]
%   rad_input_frac = p_RAD/p_input
%                  = p_RAD/(p_OHM + p_LH + p_ICRF + p_ECRH + p_NBI)
%                  = ratio of radiated power to total input power
%   rad_loss_frac  = p_RAD/p_loss
%                  = p_RAD/(p_RAD + p_cond + p_conv)
%                  = p_RAD/(p_input - dWmhd/dt)
%                  = ratio of radiated power to total loss power
%
% Author: Rewriten on EAST by Wang Bo 2015-12-15
% Some modifications by Robert Granetz 2015/10/11
% Alex Tinguely 15-09-16
% 2016/04/07  R. Granetz - Added rad_loss_frac.  To calculate this, we have
%                          to read in Wmhd from the EFIT tree
% 2017/12/11  R. Granetz - Made some changes to the ECRH section.  The ECRH
%                          power will now be read from the ANALYSIS tree.
%                          (Note that the ECRH data in the analysis tree is
%                          exactly the same as in the ECRH_EAST tree.)
%                          Also, baseline offsets will be subtracted.  And
%                          the ECRH power will now be written into the
%                          disruption warning database, as well as p_nbi.

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

% ------------------------------------
% Get lower hybrid power
%
% NOTE: the timebase for the LH power signal does not extend over the full
% time span of the discharge.  Therefore, when interpolating the LH power
% signal onto the "timebase" array, the LH signal has to be extrapolated
% with zero values.  This is an option in the 'interp1' routine.  If the
% extrapolation is not done, then the 'interp1' routine will assign NaN
% (Not-a-Number) values for times outside the LH timebase, and the NaN's
% will propagate into p_input, rad_input_frac, and rad_loss_frac, which is
% not desirable.

[shotopened, status] = mdsopen('east_1', double(shot));

if (mod(status,2)==0);
  p_LH = zeros(size(timebase));
else;
  [lh_time1, status] = mdsvalue('dim_of(\plhi1)'); % [s]
  if (mod(status,2)==0 || length(lh_time1) == 1);
    p_lh1 = zeros(size(timebase));
  else;  
    [p_lh1, status] = mdsvalue('\plhi1'); %   LHW 2.45GHz data, unit: KW 
    if (mod(status,2)==0);
      p_lh1 = zeros(size(timebase));
    else;  
      p_lh1 = 1.e3 * mdsvalue('\plhi1'); %   LHW 2.45GHz dat, convert KW to W  
      p_lh1 = interp1(lh_time1, p_lh1, timebase_column, 'linear', 0);
      if (size(timebase,2) > 1);  % If "timebase" is a row vector, then make
        p_lh1 = transpose(p_lh1); % p_lh1 into a row vector.
      end; 
    end; 
  end;
  
  [lh_time2, status] = mdsvalue('dim_of(\plhi2)'); % [s]
  if (mod(status,2)==0 || length(lh_time2) == 1);
    p_lh2 = zeros(size(timebase));
  else;  
    [p_lh2, status] = mdsvalue('\plhi2'); %   LHW 4.6GHz data, unit: KW 
    if (mod(status,2)==0);
      p_lh2 = zeros(size(timebase));
    else;
      p_lh2 = 1.e3 * mdsvalue('\plhi2'); %   LHW 4.6GHz data, convert KW to W
      p_lh2 = interp1(lh_time2, p_lh2, timebase_column, 'linear', 0);
      if (size(timebase,2) > 1);  % If "timebase" is a row vector, then make
        p_lh2 = transpose(p_lh2); % p_lh2 into a row vector.
      end; 
    end;
  end;
  
  p_LH = p_lh1 + p_lh2; % Units: W
       
  mdsclose;  
end;  
  
% ------------------------------------
% Get ECRH power

[shotopened, status] = mdsopen('analysis', double(shot));

if (mod(status,2)==0);
  p_ECRH = zeros(size(timebase));
else;
  [p_ecrh, status] = mdsvalue('\pecrh1i'); % unit [KW]
  if (mod(status,2)==0);
    p_ECRH = zeros(size(timebase));
  else;
    p_ecrh= 1.e3 * mdsvalue('\pecrh1i'); % convert [KW] to [W]
    [ecrh_time, status] = mdsvalue('dim_of(\pecrh1i)'); % [s]
    if (mod(status,2)==0 || length(ecrh_time) == 1);
      p_ECRH = zeros(size(timebase));
    else;
      baseline_indices = find(ecrh_time < 0);
      if length(baseline_indices) >= 1;
        p_ecrh_baseline = mean(p_ecrh(baseline_indices));
        p_ecrh = p_ecrh - p_ecrh_baseline;
      end;
      p_ECRH = interp1(ecrh_time, p_ecrh, timebase_column, 'linear', 0);
      if (size(timebase,2) > 1);    % If "timebase" is a row vector, then make
        p_ECRH = transpose(p_ECRH); % p_ECRH into a row vector.
      end;
    end;
  end;
  
  mdsclose;
end;

% ------------------------------------
% Get ICRF power

[shotopened, status] = mdsopen('ICRF_EAST', double(shot));

if (mod(status,2)==0);
  p_ICRF = zeros(size(timebase));
else;
  [icrf_time1, status] = mdsvalue('dim_of(\picrfii)'); % [s]
  if (mod(status,2)==0 || length(icrf_time1)==1);
    p_icrfii = zeros(size(timebase));
  else;
    [p_icrfii, status] = mdsvalue('\picrfii'); % unit: KW 
    if (mod(status,2)==0);
      p_icrfii = zeros(size(timebase));
    else;   
      p_icrfii = 1.e3 * mdsvalue('\picrfii'); % convert KW to W
      p_icrfii = interp1(icrf_time1, p_icrfii, timebase_column, 'linear', 0); 
      if (size(timebase,2) > 1);    % If "timebase" is a row vector, then make
        p_icrfii = transpose(p_icrfii); % p_icrfii into a row vector.
      end;
    end;
  end;
   
  [icrf_time2, status] = mdsvalue('dim_of(\picrfir)'); % [s]
  if (mod(status,2)==0 || length(icrf_time2)==1);
    p_icrfir = zeros(size(timebase));
  else;
    [p_icrfir, status] = mdsvalue('\picrfir'); % unit: KW
    if (mod(status,2)==0);
      p_icrfir = zeros(size(timebase));
    else; 
      p_icrfir = 1.e3 * mdsvalue('\picrfir'); % convert KW to W
      p_icrfir = interp1(icrf_time2, p_icrfir, timebase_column, 'linear', 0); 
      if (size(timebase,2) > 1);    % If "timebase" is a row vector, then make
        p_icrfir = transpose(p_icrfir); % p_icrfir into a row vector.
      end; 
    end;
  end;
     
  [icrf_time3, status] = mdsvalue('dim_of(\picrfbi)'); % [s]
  if (mod(status,2)==0 || length(icrf_time3)==1);
    p_icrfbi = zeros(size(timebase));
  else;
    [p_icrfbi, status] = mdsvalue('\picrfbi'); % unit: KW 
    if (mod(status,2)==0);
      p_icrfbi = zeros(size(timebase));
    else;   
      p_icrfbi = 1.e3 * mdsvalue('\picrfbi'); % convert KW to W
      p_icrfbi = interp1(icrf_time3, p_icrfbi, timebase_column, 'linear', 0); 
      if (size(timebase,2) > 1);    % If "timebase" is a row vector, then make
        p_icrfbi = transpose(p_icrfbi); % p_icrfbi into a row vector.
      end;
    end;
  end;
     
  [icrf_time4, status] = mdsvalue('dim_of(\picrfbr)'); % [s]
  if (mod(status,2)==0 || length(icrf_time4)==1);
    p_icrfbr = zeros(size(timebase));
  else;
    [p_icrfbr, status] = mdsvalue('\picrfbr'); % unit: KW 
    if (mod(status,2)==0);
      p_icrfbr = zeros(size(timebase));
    else;   
      p_icrfbr = 1.e3 * mdsvalue('\picrfbr'); % convert KW to W
      p_icrfbr = interp1(icrf_time4, p_icrfbr, timebase_column, 'linear', 0); 
      if (size(timebase,2) > 1);    % If "timebase" is a row vector, then make
        p_icrfbr = transpose(p_icrfbr); % p_icrfbr into a row vector.
      end;
    end;
  end;
   
  p_ICRF = (p_icrfii - p_icrfir) + (p_icrfbi - p_icrfbr); % [W]
  
  mdsclose;
end;

% ------------------------------------
% Get NBI power

[shotopened, status] = mdsopen('nbi_east', double(shot));

if (mod(status,2)==0);
  p_NBI = zeros(size(timebase));
else;
  [nbi_time1, status] = mdsvalue('dim_of(\pnbi1lsource)'); % [s]
  if (mod(status,2)==0 || length(nbi_time1) == 1);
    p_nbi1l = zeros(size(timebase));
  else;
    [p_nbi1l,status] = mdsvalue('\pnbi1lsource'); % unit [KW] 
    if (mod(status,2)==0);
      p_nbi1l = zeros(size(timebase));
    else; 
      p_nbi1l = 1.e3 * mdsvalue('\pnbi1lsource'); % convert [KW] to [W]
      p_nbi1l = interp1(nbi_time1, p_nbi1l, timebase_column, 'linear', 0);
      if (size(timebase,2) > 1);    % If "timebase" is a row vector, then make
        p_nbi1l = transpose(p_nbi1l); % p_nbi1l into a row vector.
      end;
    end;
  end;
   
  [nbi_time2, status] = mdsvalue('dim_of(\pnbi1rsource)'); % [s]
  if (mod(status,2)==0 || length(nbi_time2) == 1);
    p_nbi1r = zeros(size(timebase));
  else;
    [p_nbi1r,status] = mdsvalue('\pnbi1rsource'); % unit [KW]
    if (mod(status,2)==0);
      p_nbi1r = zeros(size(timebase));
    else; 
      p_nbi1r = 1.e3 * mdsvalue('\pnbi1rsource'); % convert [KW] to [W]
      p_nbi1r = interp1(nbi_time2, p_nbi1r, timebase_column, 'linear', 0);
      if (size(timebase,2) > 1);    % If "timebase" is a row vector, then make
        p_nbi1r = transpose(p_nbi1r); % p_nbi1r into a row vector.
      end;
    end;
  end;
    
  [nbi_time3, status] = mdsvalue('dim_of(\pnbi2lsource)'); % [s]
  if (mod(status,2)==0 || length(nbi_time3) == 1);
    p_nbi2l = zeros(size(timebase));
  else;
    [p_nbi2l,status] = mdsvalue('\pnbi2lsource'); % unit [KW]
    if (mod(status,2)==0);
      p_nbi2l = zeros(size(timebase));
    else; 
      p_nbi2l = 1.e3* mdsvalue('\pnbi2lsource'); % convert [KW] to [W]
      p_nbi2l = interp1(nbi_time3, p_nbi2l, timebase_column, 'linear', 0);
      if (size(timebase,2) > 1);    % If "timebase" is a row vector, then make
        p_nbi2l = transpose(p_nbi2l); % p_nbi2l into a row vector.
      end;
    end;
  end;
    
  [nbi_time4, status] = mdsvalue('dim_of(\pnbi2rsource)'); % [s]
  if (mod(status,2)==0 || length(nbi_time4) == 1);
    p_nbi2r = zeros(size(timebase));
  else;
    [p_nbi2r,status] = mdsvalue('\pnbi2rsource'); % unit [KW]
    if (mod(status,2)==0);
      p_nbi2r = zeros(size(timebase));
    else; 
      p_nbi2r = 1.e3 * mdsvalue('\pnbi2rsource'); % convert [KW] to [W]
      p_nbi2r = interp1(nbi_time4, p_nbi2r, timebase_column, 'linear', 0);
      if (size(timebase,2) > 1);    % If "timebase" is a row vector, then make
        p_nbi2r = transpose(p_nbi2r); % p_nbi2r into a row vector.
      end;
    end;
  end;

  p_NBI = p_nbi1l + p_nbi1r + p_nbi2l + p_nbi2r; % [W]
    
  mdsclose;
end;


% ------------------------------------
% Get ohmic power

 p_OHM = get_P_ohm_EAST(shot, timebase); %[W]

% ------------------------------------
% Radiated power

[shotopened, status] = mdsopen('prad_east', double(shot));
if (mod(status,2)==0);
  p_RAD = NaN(size(timebase));
else;
% [p_rad, status]  = mdsvalue('\pradtot_AXUV'); % [W]
  [p_rad, rad_time, status] = Prad_bulk_xuv2014_2016(shot); % original from
                                                            % Duan Yanmin, and
                                                            % modified by RSG
  if (mod(status,2)==0);
    p_RAD = NaN(size(timebase));
  else; 
    p_RAD = interp1(rad_time, p_rad, timebase_column, 'linear', 0);  
    if (size(timebase,2) > 1);  % If "timebase" is a row vector, then make
      p_RAD = transpose(p_RAD); % p_RAD into a row vector.
    end;
  end;
  mdsclose;
end;

% Get Wmhd from EFIT, calculate dWmhd/dt, and interpolate onto 'timebase'.
% This signal is needed to calculate 'rad_loss_frac'.

[shotopened, status] = mdsopen('efit18', double(shot));
if (mod(status,2) == 0);
  dWmhd_dt = NaN(size(timebase));
else;
  [Wmhd, status] = mdsvalue('\efit_aeqdsk:wplasm');
  if (mod(status,2) == 0 || length(Wmhd) <= 1);
    dWmhd_dt = NaN(size(timebase));
  else;
    efittime = mdsvalue('dim_of(\efit_aeqdsk:wplasm)');
    [efittime, indx, ~] = unique(efittime);
    Wmhd = Wmhd(indx);
    dWmhd_dt = gradient(Wmhd, efittime);
    dWmhd_dt = interp1(efittime, dWmhd_dt, timebase_column, 'linear');  
    if (size(timebase,2) > 1);        % If "timebase" is a row vector, then
      dWmhd_dt = transpose(dWmhd_dt); % make dWmhd_dt into a row vector.
    end;
  end;
end;
mdsclose;

% ------------------------------------

p_input = p_OHM + p_LH + p_ICRF + p_ECRH + p_NBI; % Total input power [W]
rad_input_frac = p_RAD ./ p_input;
rad_loss_frac  = p_RAD ./ (p_input - dWmhd_dt);

end
