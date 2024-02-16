function [duration, Ip_max] = end_of_current(Ip, Ip_time, threshold);

% Returns the time that a shot terminates on the DIII-D tokamak by
% examining the Ip signal.
%
% The optional input parameter, "threshold", is the minimum abs(current)
% required to be considered non-zero.  If not specified, it defaults
% to 100 kA
%
% Revision/history:
%  2015/01/28 - RSG; Created; based on C-Mod and EAST routines

if (~exist('threshold','var')); threshold=1.e5; end;
if (threshold < 0); threshold=abs(threshold); end; 

% Subtract baseline offset

baseindices=find(Ip_time <= 0.0);
if (numel(baseindices) > 0);
  baseline=sum(Ip(baseindices))/length(baseindices);
else;
  baseline=0;
end;
IP = Ip - baseline;

% First determine if there was any finite plasma current on this shot.  If
% not, then the shot was a "no plasma" shot, and the end-of-shot is set
% to 0 s.  Note that I do not assume uniform time sampling.

finite_indices = find((Ip_time >= 0.0) & (abs(IP) > threshold));
if (length(finite_indices) == 0);
  duration = 0;
  Ip_max = 0;
  return;
else;
  dt = diff(Ip_time);
  duration = sum(dt(finite_indices(1:end-1)));
  if (duration < 0.1); % Assumes < 100 ms is not a bona fide plasma
   duration = 0;
   Ip_max = 0;
   return;
  end;
end;

% Now determine the polarity of the plasma current.  I don't want to just
% blindly use abs(Ip), since polarity reversals may exist due to faulty
% instrumentation, and taking the absolute value would not be the same as
% simple negation.  I don't assume uniform time sampling either.

polarity = sign(trapz(Ip_time(finite_indices), IP(finite_indices)));

Ip_upright = IP * polarity;

% Find all the times that Ip is greater than the threshold.  The largest
% time value is the end of current.  But also check to see if the plasma
% current signal has been digitized for long enough to capture the end of
% the discharge.  If not, negate the end of shot value to indicate that the
% actual value cannot be determined, but it is longer than abs(value).

indices = find((Ip_upright >= threshold) & (Ip_time > 0));
duration = Ip_time(max(indices));

if (max(indices) == numel(Ip_time)); % Ip still finite at end of sampling
  duration = -duration;
end;

% Find Ip_max (with correct polarity)

Ip_max = max(Ip_upright) * polarity;

% All done
