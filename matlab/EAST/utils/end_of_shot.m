function [eos, Ip_max] = end_of_shot(shot,threshold);

% Returns the time that a shot terminates on the EAST tokamak by examining
% the Ip signal.  If the MDSplus tree containing the plasma current signal
% for the desired shot is already open when this routine is called, no
% additional mdsopen or mdsclose commands are issued, so the tree will
% remain open after the routine exits.  If no MDSplus tree is open, then
% this routine will do an mdsopen, and then an mdsclose, so the tree will
% NOT remain open upon exiting.
%
% The optional input parameter, "threshold", is the minimum abs(current)
% required to be considered non-zero.  If not specified, it defaults
% to 100 kA
%
% Revision/history:
%  2013/12/27 - RSG; Created
%  2014/01/07 - RSG; Added internal mdsopen & mdsclose procedures, and
%                      added the 'shot' parameter 
%                    Added requirement that duration > 0.1 s, where
%                      'duration' is the total integrated time for which
%                      Ip > threshold.  This eliminates a number of shots
%                      having bizarre waveforms, such as # 44196.
%  2014/01/08 - RSG; Added 'Ip_max' optional output parameter

eos = NaN;    % This default value is returned if the Ip data is unavailable
Ip_max = NaN; % This default value is returned if the Ip data is unavailable
if (~exist('threshold','var')); threshold=1.e5; end;
if (threshold < 0); threshold=abs(threshold); end; 

% Check if MDS tree is already open.  If not, open it.

mdsconnect('mds.ipp.ac.cn'); % no penalty for unnecessary invocation
[treename_open, status] = mdsvalue('getdbi("name")');
[treeshot_open, status] = mdsvalue('getdbi("shot")');
if (mod(status,2)~=0);
  if (~strcmpi(treename_open,'PCS_EAST') || (treeshot_open ~= shot));
    [shot_opened, status] = mdsopen('pcs_east',double(shot));
    if (mod(status,2)==0);
      return;
    end;
  end;
else;
  [shot_opened, status] = mdsopen('pcs_east',double(shot));
  if (mod(status,2)==0);
    return;
  end;
end;

% Read in plasma current

[Ip_time, mds_status_timebase]=mdsvalue('dim_of(\PCRL01)');
[Ip,      mds_status_signal  ]=mdsvalue(       '\PCRL01' );
if ((mod(mds_status_timebase,2)==0) || (mod(mds_status_signal,2)==0));
  fprintf(1,'Error reading Ip data for specified shot\n\n');
  if (exist('shot_opened')); mdsclose; end;
  return;
end;

% Subtract any baseline offset

baseindices=find(Ip_time <= -5.8); % time before any power supplies turn on
if (numel(baseindices) > 0);
  baseline=sum(Ip(baseindices))/length(baseindices);
  Ip = Ip - baseline;
end;

% First determine if there was any finite plasma current on this shot.  If
% not, then the shot was a "no plasma" shot, and the end-of-shot is set
% to 0 s.

finite_indices = find((Ip_time >= 0.0) & (abs(Ip) > threshold));
if (length(finite_indices) == 0);
  eos = 0;
  Ip_max = 0;
  if (exist('shot_opened')); mdsclose; end;
  return;
else;
  dt = diff(Ip_time);
  duration = sum(dt(finite_indices(1:end-1)));
  if (duration < 0.1);
   eos = 0;
   Ip_max = 0;
   if (exist('shot_opened')); mdsclose; end;
   return;
  end;
end;

% Now determine the polarity of the plasma current.  I don't want to just
% blindly use abs(Ip), since polarity reversals may exist due to faulty
% instrumentation, and taking the absolute value would not be the same as
% simple negation.  I don't assume uniform time sampling either.

polarity = sign(trapz(Ip_time(finite_indices), Ip(finite_indices)));

Ip_upright = Ip * polarity;

% Find all the times that Ip is greater than the threshold.  The largest
% time value is the end-of-shot.  But also check to see if the plasma
% current signal has been digitized for long enough to capture the end of
% the discharge.  If not, negate the end-of-shot value to indicate that
% that actual value cannot be determined, but it is longer than
% abs(value).

indices = find((Ip_upright >= threshold) & (Ip_time > 0));
eos = Ip_time(max(indices));

if (max(indices) == numel(Ip_time)); % Ip still finite at end of sampling
  eos = -eos;
end;

% Find Ip_max (with correct polarity)

Ip_max = max(Ip_upright) * polarity;

% All done.  Close the tree if necessary, and return the values.

if (exist('shot_opened')); mdsclose; end;
