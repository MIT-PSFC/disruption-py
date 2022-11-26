function [status, shotopened, disrupt_flag, t_disrupt, Ip0, max_dIdt] = ...
  test_for_disruption(shot);

% This routine determines whether or not a given D3D plasma shot
% terminated in a disruption.  If it finds that the shot disrupted, it also
% returns the time of the disruption (t_disrupt), the plasma current just
% prior to the disruption (Ip0), and the maximum dI/dt during the current
% quench (max_dIdt).
% 
% The following criteria (partially based on those used by Robert Granetz
% and Steve Wolfe for Alcator C-Mod) are used to determine disruptivity.
% The parameters have been adjusted for the D3D tokamak.
%
% 1) shot duration > 0.5 s (reject very short shots)
% 2) abs(Ip_max) > 0.1 MA (reject very low current shots)
% 3) Ip0/Ip_max > 0.33  (reject shots that disrupt late in current rampdown)
% 4) Ip0/max_dIdt < 0.050 s (reject shots with relatively slow current decay)
% 5) |Ip_final| < 100 kA within 0.15 s of t_disrupt (reject minor disruptions)
% 6) abs(Ip0) > 0.1 MA (reject very low current disruptions)
%
% Revision/history:
%  2015/01/26 - RSG; Created.

% Set values for test parameters (tuned for D3D)

duration_min = 0.5; % minimum shot duration (seconds)
Ip_threshold = 0.1e6; % minimum plasma current (amperes)
tau_CQ_max = 0.050; % upper limit on current quench (CQ) decay time (seconds)
Ip_final_max = 100e3; % upper limit on post-disruption current

% Set default return values

disrupt_flag=NaN;
t_disrupt=NaN;
Ip0=NaN;
max_dIdt=NaN;

mdsconnect('atlas.gat.com');

[shotopened, mds_status]=mdsopen('d3d',shot); % Make sure the shot tree
if (mod(mds_status,2)==0);                    % exists.
  shotopened = shot;
  fprintf(1,'  Unable to open D3D tree for shot%7i\n', shotopened);
  status = 0;
  return;
else; % Read in plasma current.
  [Ip_time, mds_status_timebase] = ...
    mdsvalue(['dim_of(ptdata("ip",' num2str(shot) '))']); % milliseconds
  [Ip, mds_status_signal  ] = ...
    mdsvalue(['ptdata("ip",' num2str(shot) ')']); % amperes
  mdsclose;
  if ((mod(mds_status_timebase,2)==0) || (mod(mds_status_signal,2)==0));
    fprintf(1,'  Error reading Ip data for shot%7i\n', shotopened);
    status = 0;
    return;
  else;
    status = 1;
    Ip_time = Ip_time/1.e3; % Convert to seconds
  end;
end;

% Subtract baseline offset

baseindices=find(Ip_time <= 0.0);
if (numel(baseindices) > 0);
  baseline=sum(Ip(baseindices))/length(baseindices);
  Ip = Ip - baseline;
end;

% Determine duration of shot by calling the routine "end_of_current".

duration = end_of_current(Ip, Ip_time, Ip_threshold);

% Test candidate shot against criteria #1 and #2.  (Criterion #2 is
% implicitly handled by the "end_of_shot" routine, which returns zero for
% the duration if Ip < Ip_threshold.)

if (duration < duration_min);
  disrupt_flag = 0;
  return;
end;

% Determine the polarity of the plasma current.  I don't want to just
% blindly use abs(Ip), since polarity reversals may exist due to faulty
% instrumentation, and taking the absolute value would not be the same as
% simple negation.

time_indices = find((Ip_time > 0) & (Ip_time < duration));
polarity = sign(trapz(Ip_time(time_indices), Ip(time_indices)));

Ip_upright = Ip * polarity;

% Find maximum plasma current, but try to avoid any disruption spike (if
% the shot disrupted; we still don't know that yet.)

time_indices = find((Ip_time > 0) & (Ip_time < duration - 0.050));
Ip_max = max(Ip_upright(time_indices)) * polarity;

% Find plasma current during the period just before the end of the
% discharge.  This will be the pre-disruption plasma current, if the shot
% actually disrupted.

time_indices = find((Ip_time > duration - 0.06) & ...
                    (Ip_time < duration - 0.04));

% Shot 124948 has a screwed up timebase.  The following statement deals
% generally with that special case.

if (numel(time_indices)==0);
  status = 0;
  fprintf(1,'  Timebase for shot%7i is invalid\n', shotopened);
  return;
end;

candidate_Ip0 = mean(Ip_upright(time_indices)) * polarity;

% Apply criterion #3 to reject shots that disrupt late in the normal
% current rampdown

if (candidate_Ip0 / Ip_max < 0.33);
  disrupt_flag=0;
  return;
end;

% Now compute dI/dt during the latter part of the discharge.  If it's
% faster than the specified limit, then this shot is a candidate for being
% declared a disruption, and the time of maximum abs(dI/dt) is the
% disruption time.

time_indices = find((Ip_time > duration - 0.05) & ...
                    (Ip_time < duration + 0.05));

dI_upright = diff(Ip_upright(time_indices));
dt = diff(Ip_time(time_indices));
dIdt_upright = dI_upright./dt;
[candidate_max_dIdt, indx] = min(dIdt_upright);
candidate_max_dIdt = candidate_max_dIdt * polarity;
candidate_t_disrupt = Ip_time(time_indices(indx));

% Now test the candidate against criterion #4

if (-candidate_Ip0/candidate_max_dIdt > tau_CQ_max);
  disrupt_flag=0;
  return;
end;

% Now test the candidate against criterion #5

time_indices = find((Ip_time > candidate_t_disrupt) & ...
                    (Ip_time < candidate_t_disrupt + 0.15));
Ip_final = abs(min(Ip_upright(time_indices)));
if (Ip_final > Ip_final_max);
  disrupt_flag=0;
  return;
end;

% Now test the candidate against criterion #6

if (abs(candidate_Ip0) < Ip_threshold);
  disrupt_flag=0;
  return;
end;

% Okay, we have a disruption!

disrupt_flag = 1;
t_disrupt = candidate_t_disrupt;
Ip0 = candidate_Ip0;
max_dIdt = candidate_max_dIdt;
