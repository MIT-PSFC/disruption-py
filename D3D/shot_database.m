if (~exist('EAST_shot_list.txt','file'));
  fileid=fopen('EAST_shot_list.txt','wt'); % create ASCII datafile
  fprintf(fileid, ...
    ['   Shot         Ip0       end of shot          timestamp     \n' ...
     '                (A)           (s)                            \n' ...
     ' -------     --------     ----------      -------------------\n']);
  fclose(fileid);
end;

Ip_threshold = 0.1e6; % threshold for finite plasma current is 0.1 MA

if (~exist('startshot'));
  startshot = input('Enter starting shot number: ');
  if (isempty(startshot));
    clear startshot;
    return;
  end;
end;
if (~exist('endshot'));
  endshot = input('Enter ending shot number: ');
  if (isempty(endshot));
    endshot = startshot + 1000;
  end;
end;
for shot = startshot : endshot;
  fprintf(1,'Processing shot %7i \n', shot);
  [eos, Ip0] = end_of_shot(shot, Ip_threshold);
%  Get data and time of shot
  timestamp = NaN;
  [dummy, mds_status] = mdsopen('east',shot);
  if (mod(mds_status,2)==1);
    [createtime, mds_status_signal] = mdsvalue('\ipm:createtime');
    if (mod(mds_status_signal,2)==1);
      date_array = datevec(createtime);        % Convert date/time format
      year   = num2str(date_array(1),'%4i');   % into SQL-compatible format.
      month  = num2str(date_array(2),'%02i');
      day    = num2str(date_array(3),'%02i');
      hour   = num2str(date_array(4),'%02i');
      minute = num2str(date_array(5),'%02i');
      second = num2str(date_array(6),'%02i');
      timestamp = [year '-' month '-' day ' ' hour ':' minute ':' second];
     end;
    mdsclose;
  end;
  fileid=fopen('EAST_shot_list.txt','at'); % open file w/append access
  fprintf(fileid,' %7i     %8i     %10.4f      %19s\n', ...
                  shot, round(Ip0),  eos,  timestamp);
  fclose(fileid);
end;
