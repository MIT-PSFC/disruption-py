function [n_equal_1_mode, n_equal_1_normalized] = ...
  get_n1_bradial_d3d(shot, timebase);

% This routine obtaines the n=1 Bradial information from the ESLD coils (units of Tesla).
% It also calculates the normalized Bradial, dividing by the toroidal B-field.
% For shots after 176030, a manual recalculation of the Bradial was needed;
% B. Sammuli reprocessed 880 shots that were missing the hardware compensations
% for the F-coils.[176030-176912], first good shot with calibrated Bradial is
% 176913.

% Inputs:
%   shot = shot number
%   timebase = array of desired time values
%
% Outputs:
%   n_equal_1_mode = bradial n1 component calculated by ESLD coils [T]
%   n_equal_1_normalized = n_equal_1_mode / Btor
%
%
% Author: Cristina Rea and Robert Granetz   25 Jul 2018
% Revision/history:
%  2018/08/16 - RSG: inserted code to check if bradial data exists in the
%                    NetCDF data for the requested shot, if that shot is
%                    within the range that's supposed to be in the NetCDF
%                    dataset.  If there are no data for the requested shot,
%                    then return NaN's. 
%  2019/04/05 - CR:  This routine will not find any dusbradial or onsbradial
%                    data in the "standard" (atlas) tree paths for certain
%                    discharges. In September 2017, I went through the 
%                    existing 2014-2017 SQL database discharges and found 
%                    that 3814 shots were missing the bradial calculation.
%                    B. Sammuli recomputed the missing bradial from the 
%                    existing ESLDs and stored the info as MDS data, in 
%                    tree structures for the specified shot numbers.
%                    The actual data is stored in 
%                    /fusion/projects/disruption_warning/recalc_bradial
%                    and /fusion/projects/codes/pcs/data/mdsdata/bradial
%FIXME
%                    Implement check that if shot number is in the list:
%                    disruption_warning/recalc_bradial/shotlist_dusbradial_reprocessed.txt
%                    then the correct data source for the bradial is in
%                    disruption_warning/recalc_bradial/bradial

% The input array, "timebase", can be either a row vector (1 x n) or a
% column vector (n x 1).  We desire to have the output arrays match the
% shape of "timebase".  However, Matlab calls to the routine "mdsvalue"
% return column vectors for 1-D signals, so it is simpler to work with
% column vectors within this routine, and then, if necessary, convert the
% outputs to row vectors just before exiting this routine.  So the first
% step is to create a copy of "timebase" that is guaranteed to be a column
% vector.

if (size(timebase,1) > 1);
  timebase_column = timebase;
else;
  timebase_column = transpose(timebase);
end;

% Initialize all output arrays to NaN (Not-a-Number) column vectors

n_equal_1_mode = NaN(size(timebase));
n_equal_1_normalized = NaN(size(timebase));

mdsconnect('atlas.gat.com');
% Get data from NetCDF file if shot is in [176030-176912]
% otherwise get bradial from DUD or ONFR systems

if shot >= 176030 & shot <= 176912; % First shot with missing compensation for bradial calculations
  fprintf(1,['Reading bradial from NetCDF for shot%7i\n'], shot);
  filename = '/fusion/projects/disruption_warning/matlab_programs/recalc.nc';
  ncid = netcdf.open(filename,'NOWRITE');% open NetCDF
  brad = ncread(filename,'dusbradial_calculated');
  n1time = ncread(filename,'times'); %[ms]
  n1time = n1time*1.e-3; %[s] 
  shots = ncread(filename,'shot');
  indx = find(shots == shot);              % Modified by R. Granetz
  if length(indx) == 1;                    % on 2018/08/16 to check
    dusbradial = brad(:,indx)*1.e-4;%[T]   % for, and deal with,
  else;                                    % shots that are missing
    dusbradial = NaN(size(brad,1), 1);     % in the NetCDF dataset.
  end;
  netcdf.close(ncid);
else; % if statement to check if shot is in shotlist_dusbradial_reprocessed.txt
  fprintf(1,['Reading bradial from MDS+ for shot%7i\n'], shot);
  [n1time,status] = mdsvalue(['dim_of(ptdata("dusbradial", ' num2str(shot) '))']);%[ms]
  if (mod(status,2) == 1 && length(n1time) >= 4);%1
    fprintf(1,'dusbradial\n')
    n1time = n1time/1.e3; % convert ms to s
    dusbradial = mdsvalue(['ptdata("dusbradial", ' num2str(shot) ')'])*1.e-4; %[T]
    mdsclose;
  else;
  % If there's no data in dusbradial (DUD system)
  % check if ONFR was on and calculating the bradial
    [n1time,status] = mdsvalue(['dim_of(ptdata("onsbradial", ' num2str(shot) '))']); 
    if (mod(status,2) == 1 && length(n1time) >= 4);%2
      fprintf(1,'onsbradial\n')
      n1time = n1time/1.e3; % convert ms to s
      dusbradial = mdsvalue(['ptdata("onsbradial", ' num2str(shot) ')'])*1.e-4; %[T]
      mdsclose;
    else; % return NaNs if no bradial was calculated
      if (size(timebase,2) > 1);
      n_equal_1_mode = transpose(n_equal_1_mode);
      n_equal_1_normalized = transpose(n_equal_1_normalized);
      end;
      return;
    end;%2
  end;%1
end;%0
n_equal_1_mode = interp1(n1time, dusbradial, timebase_column, 'linear'); % [T]

% Next, get the toroidal field Btor, which is needed to calculate the 
% n_equal_1_normalized

[Btime, status] = mdsvalue(['dim_of(ptdata("bt", ' num2str(shot) '))']);% ms
if (mod(status,2) == 1);
  Btime = Btime/1.e3; % convert ms to s
  Btor = mdsvalue(['ptdata("bt", ' num2str(shot) ')']); %[T]
  mdsclose;
  Btor = interp1(Btime, Btor, timebase_column, 'linear');
  n_equal_1_normalized = n_equal_1_mode ./ abs(Btor);
else;
  if (size(timebase,2) > 1);
    n_equal_1_mode = transpose(n_equal_1_mode);
    n_equal_1_normalized = transpose(n_equal_1_normalized);
  end;
  return;
end;

end
