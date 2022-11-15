function disruption_database(shotlist, varargin);

% 'shotlist' is a Matlab array of shot numbers.  It can be a single shot,
%   or a sequence of shots [159000:160000], or a list of specific shots
%   [159490, 159499, 159524], or combinations of these [159495, 159524,
%   159600:160000].  Shot #0 refers to the current DIII-D shot.
%
% 'varargin' allows for explicit specification of the fields that will be
%   updated in the disruption table.  The varargin parameter(s) must be
%   character strings.  If varargin is not specified, then all fields will
%   be updated.
%
% Examples of usage:
%
%  disruption_database(0) % Analyze current DIII-D shot.  If it disrupted,
%                         % insert data into all fields (a new record).
%
%  disruption_database(159000:162000)
%                         % Analyze all DIII-D shots from #159000 to 162000.
%                         % For each shot that disrupted, insert data into
%                         % all fields.  If the shot already exists in the
%                         % database, then insert data into all fields
%                         % except the shot number.  (It's the primary key.)
%
%  disruption_database(159000:162000, 'i_halo_l', 'i_halo_u', 'tpf')
%                         % Analyze all DIII-D shots from #159000 to 162000.
%                         % For each shot that disrupted, insert data into
%                         % fields 'i_halo_u', 'i_halo_l', and 'tpf' only.
%
% Revision/history:
%  2015/01    - RSG; Started (based on my program for EAST)
%  2015/05    - RSG; Get program running by using primitive JDBC driver calls
%                    because DIII-D does not have the Matlab Database Toolkit

fieldlist = lower(strtrim(varargin)); % lowercase; remove spaces

% Since Matlab at DIII-D does not have the Database Toolkit, the higher
% level routines to connect, query, and insert/update data records into the
% database are not available.  Therefore, low-level calls to the JDBC
% driver need to be used instead.  What a pain ... and totally illogical,
% unreadable, inefficient, etc.

% Connect to the D3DRDB sql database and determine how many fields are in
% the disruption table.

db = set_database('d3drdb'); % My routine to connect to database

%result = fetch(db, 'describe disruptions');
%fields_in_db = result(:,1);

meta = db.getMetaData();
rs = meta.getColumns( {}, {}, 'disruptions', {} );
fields_in_db={};
while (rs.next());
  fields_in_db = [fields_in_db, rs.getObject(4)];
end;

if (length(fieldlist) == 0);
  nfields = length(fields_in_db);  
else;
  nfields = length(fieldlist);
  not_okay = ~ismember(fieldlist, fields_in_db);
  if (sum(not_okay) > 0);
    fprintf(1, '\nThese invalid field names will be ignored:\n');
    fprintf(1, '   %s\n', fieldlist{1, not_okay});
    fprintf(1, '\n');
  end;
end;

% Get list of shots that are already in the disruption table

%result = fetch(db, 'select shot from disruptions order by shot');
%shots_already_in_db = cell2mat(result);

% Because I have to use the low-level brain-dead calls to the JDBC driver
% (since Matlab at DIII-D does not have the Database Toolkit), I first have
% to find out how many shot records are already in the disruption database
% so that I can pre-allocate an array to hold those shot numbers.  What a
% pain!

ps = db.prepareStatement('select count(*) from disruptions');
rs = ps.executeQuery();
rs.next();
nshots_already_in_db = rs.getObject(1);
shots_already_in_db = zeros(1,nshots_already_in_db); % Allocate array

% Now that I have allocated the array, I can go ahead and fill it in with
% the shot numbers that already have records in the disruption database.

ps = db.prepareStatement('select shot from disruptions order by shot');
rs = ps.executeQuery();
for ishot = 1:nshots_already_in_db;
  rs.next();
  shots_already_in_db(ishot) = rs.getObject(1);
end;

% Okay, we're all set up.  Now loop through the list of shots.  Determine
% whether each shot disrupted.  If the shot disrupted, then gather the
% disruption-relevant data, and insert it into the disruption database.

for shot = shotlist;
  [status, shotopened, flag, tdis, ip, didt] = test_for_disruption(shot);
  if (~status);
    fprintf(1,'Shot%8i could not be evaluated for disruption\n', shotopened);
    continue;
  end;
  if (~flag);
    fprintf(1,'Shot%8i did not disrupt\n', shotopened);
    continue; % Skip to the next shot in the list
  else;
    fprintf(1,'Shot%8i disrupted\n', shotopened);
  end;

% Create two cell arrays, 'fields' and 'values'.  These will be used to
% enter data into the database for each specified shot.

  fields = cell(1,nfields);
  values = cell(1,nfields);

% Begin filling in the two arrays, 'fields' and 'values'.

  field_counter = 0;
  if (length(fieldlist) == 0);
    fields(1,1:4) = {'shot','ip','t_disrupt','didt'};
    values(1,1:4) = {shotopened, ip, tdis, didt};
    field_counter = field_counter + 4;
  else;
    if ismember('shot', fieldlist);
      field_counter = field_counter + 1;
      fields{1,field_counter} = 'shot';
      values{1,field_counter} = shotopened;
    end;
    if ismember('ip', fieldlist);
      field_counter = field_counter + 1;
      fields{1,field_counter} = 'ip';
      values{1,field_counter} = ip;
    end;
    if ismember('t_disrupt', fieldlist);
      field_counter = field_counter + 1;
      fields{1,field_counter} = 't_disrupt';
      values{1,field_counter} = tdis;
    end;
    if ismember('didt', fieldlist);
      field_counter = field_counter + 1;
      fields{1,field_counter} = 'didt';
      values{1,field_counter} = didt;
    end;
  end;

% Get values for all the other fields in the disruption database.  For each
% valid value, insert the field name and corresponding value into the
% 'fields' and 'values' cell arrays respectively.

% But first, create a string containing the shot number.  This will be
% needed over and over again in calls to the D3D 'ptdata' routine.

shotstring = num2str(shotopened);

% Start with timestamp (the date and time of the shot).  On DIII-D this is
% read from the 'summaries' table in the D3DRDB database.  This requires
% using primitive calls to the JDBC SQL driver, since DIII-D doesn't have
% the Matlab Database Toolkit.

  if (length(fieldlist) == 0 || ismember('timestamp', fieldlist));
    ps = db.prepareStatement(['select time_of_shot from summaries ' ...
          'where shot = ' shotstring]);
    rs = ps.executeQuery();
    if (rs.next()==1);
      timestamp = char(rs.getObject(1));
      if (length(timestamp > 19)); timestamp = timestamp(1:19); end;
      field_counter = field_counter + 1;
      fields{1,field_counter} = 'timestamp';
      values{1,field_counter} = timestamp;
    end;
  end;

% Next, get the toroidal magnetic field.

  if (length(fieldlist) == 0 || ismember('bt', fieldlist));
    [bt_sig, mds_status_signal] = mdsvalue(['ptdata("bt",' shotstring ')']);
    [bt_time, mds_status_timebase] = ...
       mdsvalue(['dim_of(ptdata("bt",' shotstring '))']);
    if (mod(mds_status_signal,2)==1 && mod(mds_status_timebase,2)==1);
      bt_time = bt_time/1.e3; % convert ms to s
      time_indices = find((bt_time > tdis - 0.06) & ...
                          (bt_time < tdis - 0.04));
      bt = mean(bt_sig(time_indices));
      field_counter = field_counter + 1;
      fields{1,field_counter} = 'bt';
      values{1,field_counter} = bt ;
    end;
  end;

% Next, get EFIT data.

  [dummy, efit_status] = mdsopen('efit01',shotopened);
  if (mod(efit_status,2)==1);
    [dummy, efit_status] = mdsvalue('tcl("set default \\efit_a_eqdsk")');
    [efit_time, efit_time_status] = mdsvalue('\atime'); % units: ms
    if (mod(efit_time_status,2)==1 ...  % If EFIT time is not in tree, then
        && (numel(efit_time)>1));       % EFIT did not run on this shot.  No
                                        % need to get any other EFIT data.
                                        % Or if EFIT time only has one
                                        % element, then EFIT crashed.
      efit_time = efit_time/1.e3; % convert ms to s
      efit_index = max(find(efit_time < tdis)); % Find last EFIT time
                                                % before the disruption

      if (length(fieldlist) == 0 || ismember('aminor', fieldlist));
        [aminor, mds_status_signal] = mdsvalue('\aminor'); % units: m
        if (mod(mds_status_signal,2)==1);
          field_counter = field_counter + 1;
          fields{1,field_counter} = 'aminor';
          values{1,field_counter} = aminor(efit_index);
        end;
      end;

      if (length(fieldlist) == 0 || ismember('raxis', fieldlist));
        [raxis, mds_status_signal] = mdsvalue('\r0'); % units: m
        if (mod(mds_status_signal,2)==1);
          field_counter = field_counter + 1;
          fields{1,field_counter} = 'raxis';
          values{1,field_counter} = raxis(efit_index);
        end;
      end;

      if (length(fieldlist) == 0 || ismember('zaxis', fieldlist));
        [zaxis, mds_status_signal] = mdsvalue('\z0'); % units: m
        if (mod(mds_status_signal,2)==1);
          field_counter = field_counter + 1;
          fields{1,field_counter} = 'zaxis';
          values{1,field_counter} = zaxis(efit_index);
        end;
      end;

      if (length(fieldlist) == 0 || ismember('wmhd', fieldlist));
        [wmhd, mds_status_signal] = mdsvalue('\wmhd'); % units: J
        if (mod(mds_status_signal,2)==1);
          field_counter = field_counter + 1;
          fields{1,field_counter} = 'wmhd';
          values{1,field_counter} = wmhd(efit_index);
        end;
      end;

      if (length(fieldlist) == 0 || ismember('topology', fieldlist));
        [topology, mds_status_signal] = mdsvalue('limloc');
        if ((mod(mds_status_signal,2)==1));
          field_counter = field_counter + 1;
          fields{1,field_counter} = 'topology';
          values{1,field_counter} = topology{efit_index};
        end; 
      end;

      if (length(fieldlist) == 0 || ismember('betap', fieldlist));
        [betap, mds_status_signal] = mdsvalue('\betap'); % dimensionless
        if (mod(mds_status_signal,2)==1);
          field_counter = field_counter + 1;
          fields{1,field_counter} = 'betap';
          values{1,field_counter} = betap(efit_index);
        end;
      end;

      if (length(fieldlist) == 0 || ismember('elong', fieldlist));
        [kappa, mds_status_signal] = mdsvalue('\kappa'); % dimensionless
        if (mod(mds_status_signal,2)==1);
          field_counter = field_counter + 1;
          fields{1,field_counter} = 'elong';
          values{1,field_counter} = kappa(efit_index);
        end;
      end;

      if (length(fieldlist) == 0 || ismember('li', fieldlist));
        [li, mds_status_signal] = mdsvalue('\li'); % dimensionless
        if (mod(mds_status_signal,2)==1);
          field_counter = field_counter + 1;
          fields{1,field_counter} = 'li';
          values{1,field_counter} = li(efit_index);
        end;
      end;

      if (length(fieldlist) == 0 || ismember('q95', fieldlist));
        [q95, mds_status_signal] = mdsvalue('\q95'); % dimensionless
        if (mod(mds_status_signal,2)==1);
          field_counter = field_counter + 1;
          fields{1,field_counter} = 'q95';
          values{1,field_counter} = q95(efit_index);
        end;
      end;

      if (length(fieldlist) == 0 || ismember('gap_inner', fieldlist));
        [gapin, mds_status_signal] = mdsvalue('\gapin'); % units: m
        if (mod(mds_status_signal,2)==1);
          field_counter = field_counter + 1;
          fields{1,field_counter} = 'gap_inner';
          values{1,field_counter} = gapin(efit_index);
        end;
      end;

      if (length(fieldlist) == 0 || ismember('gap_outer', fieldlist));
        [gapout, mds_status_signal] = mdsvalue('\gapout'); % units: m
        if (mod(mds_status_signal,2)==1);
          field_counter = field_counter + 1;
          fields{1,field_counter} = 'gap_outer';
          values{1,field_counter} = gapout(efit_index);
        end;
      end;

      if (length(fieldlist) == 0 || ismember('betan', fieldlist));
        [betan, mds_status_signal] = mdsvalue('\betan');
        if (mod(mds_status_signal,2)==1);
          field_counter = field_counter + 1;
          fields{1,field_counter} = 'betan';
          values{1,field_counter} = betan(efit_index);
        end;
      end;

%  Additional possible EFIT parameter: poh

    else;
      clear efit_index;
    end;
    mdsclose;
  end;

% Now get some non-EFIT data

  if (length(fieldlist) == 0 || ismember('direction', fieldlist));
    [zcentroid, mds_status_signal] = ...
      mdsvalue(['ptdata("vpsdfz",' shotstring ')']);
    [tzcentroid, mds_status_timebase] = ...
      mdsvalue(['dim_of(ptdata("vpsdfz",' shotstring '))']);
    if (mod(mds_status_signal,2)==1 && mod(mds_status_timebase,2)==1);
      zcentroid = zcentroid * 0.1; % convert to m
      tzcentroid = tzcentroid/1.e3; % convert ms to s
      [dummy,indx] = min(abs(tzcentroid - (tdis+0.003)));
      zc = zcentroid(indx);
      if (zc <= -0.02);
        direction = -1;
      elseif (zc >= +0.02);
        direction = +1;
      else;
        direction =  0;
      end;
      field_counter = field_counter + 1;
      fields{1,field_counter} = 'direction';
      values{1,field_counter} = direction;
    end;
  end

% Now determine the "flattop" parameter, which is 1 (true) if the shot
% disrupted during the current flattop, and 0 (false) otherwise.  Use the
% target Ip signal (pointname "iptipp") to figure out the flattop times.

  if (length(fieldlist) == 0 || ismember('flattop', fieldlist));
    [target_ip, mds_status_signal] = ...
      mdsvalue(['ptdata("iptipp",' shotstring ')']);
    [target_ip_time, mds_status_timebase] = ...
      mdsvalue(['dim_of(ptdata("iptipp",' shotstring '))']);
    if ((mod(mds_status_timebase,2)==1) && (mod(mds_status_signal,2)==1));
      target_ip_time = target_ip_time/1.e3; % convert ms to s
      max_target_ip = max(abs(target_ip));
      flattop_indices = find(abs(target_ip)/max_target_ip > 0.999);
      flattop_start = target_ip_time(min(flattop_indices));
      flattop_end   = target_ip_time(max(flattop_indices));
      flattop_duration = flattop_end - flattop_start;
      if (tdis >= flattop_start && tdis <= flattop_end+0.010 && ...
          flattop_duration >= 0.100);
       flattop = 1;
      else;
       flattop = 0;
      end;
      field_counter = field_counter + 1;
      fields{1,field_counter} = 'flattop';
      values{1,field_counter} = flattop;
    end;
  end;

% Now get information on nebar, Te, Prad, PNBI, PLH, PICRF, etcetera.

% For the disruption database, I want to select the value of each of these
% signals at a time just before the disruption.  If the EFIT data exists,
% then use the last EFIT time before the disruption.  If the EFIT data does
% not exist, then use 'tdis - 0.050'

  if (mod(efit_status,2)==1 && mod(efit_time_status,2)==1 && ...
      exist('efit_index') && ~isempty(efit_index));
    t_sample = efit_time(efit_index);
  else;
    t_sample = tdis - 0.050;
  end;

  if (length(fieldlist) == 0 || ismember('nebar', fieldlist));
    [dummy, mds_status] = mdsopen('bci',shotopened);
    if (mod(mds_status,2)==1);
      [nebar     , mds_status_signal  ]=mdsvalue('       \denv2 ');
      [nebar_time, mds_status_timebase]=mdsvalue('dim_of(\denv2)');
      if ((mod(mds_status_timebase,2)==1) && (mod(mds_status_signal,2)==1));
        nebar_time = nebar_time/1.e3; % convert ms to s
        [dummy, nebar_time_index] = min(abs(nebar_time - t_sample));
        field_counter = field_counter + 1;
        fields{1,field_counter} = 'nebar';
        values{1,field_counter} = nebar(nebar_time_index)/2.774*1.e6; % m-3
      end;
      mdsclose;
    end;
  end;

  if (length(fieldlist) == 0 || ismember('te0', fieldlist));
    [te0, mds_status_signal] = ...
      mdsvalue(['ptdata("ece1",' shotstring ')']);
    [te0_time, mds_status_timebase] = ...
      mdsvalue(['dim_of(ptdata("ece1",' shotstring '))']);
    if (mod(mds_status_signal,2)==1 && mod(mds_status_timebase,2)==1);
      te0_time = te0_time/1.e3; % convert ms to s
      [dummy, te0_time_index] = min(abs(te0_time - t_sample));
      field_counter = field_counter + 1;
      fields{1,field_counter} = 'te0';
      values{1,field_counter} = te0(te0_time_index)*1.e3; % eV
    end;
  end;

% Okay, we're done with getting the data.  Now prepare to enter the data
% into the disruption database.

  if (field_counter == 0); % If no data is available for any of the fields,
    continue;              % then skip to the next shot in the list
  end;

  fields = fields(1,1:field_counter); % Truncate these cell arrays to get
  values = values(1,1:field_counter); % rid of unassigned fields.

% Also get rid of any fields which have NaN for the corresponding data.

  test_fields = zeros(1,field_counter); % The 'isnan' function does not work
  for field_index = 1:field_counter;    % on entire cell arrays, but does
    if (isnan(values{1,field_index}));  % work on individual cell elements.
      test_fields(field_index) = NaN;   % So create a corresponding non-cell
    end;                                % array to perform the 'isnan' test.
  end;
  good_indices = ~isnan(test_fields);
  fields = fields(1,good_indices);
  values = values(1,good_indices);
  field_counter = sum(good_indices);

% Now enter the data for this shot into the disruption table in the database.
% If the shot already has a record in the table, then this is an update
% operation.  All fields EXCEPT 'shot' can be updated.  If a record for this
% shot does not already exist in the table, then this operation will insert
% a completely new record, including the 'shot' field.

  if ismember(shotopened, shots_already_in_db); % If a record for this shot
                                                % is already in the database,
                                                % then:
    if ismember('shot', fields); % Determine if 'shot' is in list of fields.
      indices = find(~ismember(fields, 'shot')); % If yes, then find indices
                                                 % of all fields that are
                                                 % not 'shot'.
      fields = fields(1, indices); % Remove 'shot' field.
      values = values(1, indices); % Remove shot value.
      field_counter = numel(indices);
    end;
    if (field_counter > 0);
%     update(db, 'disruptions', fields, values, ...
%       ['where shot = ' num2str(shotopened,'%i')]);

% Since Matlab at DIII-D doesn't have the Database Toolkit, I can't use
% the 'unpdate' routine.  Instead, I have to resort to primitive JDBC
% driver calls.  This also requires that I create a string array of all
% the fields to be inserted, with question marks ('?') as placeholders for
% all the values to be inserted.

    field_string = '';
    for field_index = 1 : field_counter;
      field_string = [field_string fields{field_index} ' = ?'];
      if (field_index ~= field_counter);
        field_string = [field_string ', ' ]; % Need comma between fields
      end;
    end;

    sql_string = ['update disruptions set ' field_string ...
                  ' where shot = ' num2str(shotopened,'%i')];
    ps = db.prepareStatement(sql_string);
    for field_index = 1 : field_counter;
      ps.setObject(field_index, values{field_index});
    end;
    rs = ps.executeUpdate();

    end;
  else;                                         % If this shot is not already
                                                % in the database, then:
    if ~ismember('shot', fields); % Determine if 'shot' is NOT in list of
                                  % fields.  The 'shot' field is required to
                                  % create a new record in disruptions table.
      fields = ['shot', fields]; % Add 'shot' field if not already in list.
      values = [shotopened, values]; % Add shot value if not already in list.
      field_counter = field_counter + 1;
    end;

%   fastinsert(db, 'disruptions', fields, values);

% Since Matlab at DIII-D doesn't have the Database Toolkit, I can't use
% the 'fastinsert' routine.  Instead, I have to resort to primitive JDBC
% driver calls.  This also requires that I create a string array of all
% the fields to be inserted, and another string array of question marks
% ('?') for all the values to be inserted.

    field_string = '('; % Start with an opening parenthesis
    value_string = '('; % Start with an opening parenthesis
    for field_index = 1 : field_counter;
      field_string = [field_string fields{field_index}];
      value_string = [value_string '?'];
      if (field_index ~= field_counter);
        field_string = [field_string ', ' ]; % Need comma between fields
        value_string = [value_string ','  ]; % Need comma between values
      else;
        field_string = [field_string ')' ];  % Need to close parentheses
        value_string = [value_string ')' ];  % Need to close parentheses
      end;
    end;

    sql_string = ['insert into disruptions ' field_string ' values ' ...
                  value_string];
    ps = db.prepareStatement(sql_string);
    for field_index = 1 : field_counter;
      ps.setObject(field_index, values{field_index});
    end;
    rs = ps.executeUpdate();

  end;

% Okay, this shot is done.  Go to the next shot in shotlist.

end;

db.close(); % Close database server connections
clear db;
mdsdisconnect;
