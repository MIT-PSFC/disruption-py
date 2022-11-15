cd('/fusion/projects/disruption_warning/matlab_programs/');

db = set_database('d3drdb');  % Connect to SQL database and return handle "db"

% Create a cell array containing the names of three of the columns (i.e.
% variables) in the "disruption_warning" table.

fields = {'shot', 'time', 'q95'};

% Construct three arrays of data corresponding to the aforementioned
% variables.  Let each array have five values, i.e. data sampled at
% five time slices.

shot = [172944, 172944, 172944, 172944, 172944];
time = [   0.5,    1.5,    2.5,    3.5,    4.5]; % seconds
q95  = [   5.0,    4.8,    NaN,    4.4,    4.2]; % Include one NaN as a test

% Assemble the 1-D arrays of data into a single 2-D cell array in the format
% required by Matlab's SQL "insert" function

values = cell(5, 3);  % 5 values in time for each of 3 variables
values(:, 1) = num2cell(shot);
values(:, 2) = num2cell(time);
values(:, 3) = num2cell( q95);

% Each record in our SQL database contains the data for a single time slice
% of a single shot.  The "shot" and "time" columns are primary keys,
% i.e. they are unique for each record.  Since the specified shot (172944)
% does not yet have any time slices in the existing database, we need to
% insert new records for these time slices, as opposed to simply updating
% the "q95" values in existing records.  Therefore we use Matlab's SQL
% "insert" function rather than the SQL "update" function.  (Note: with the
% java interface, the "insert" function actually calls the "fastinsert"
% function.)

for itime = 1:5;
  fprintf(1, 'Now calling Matlab''s SQL "insert" function\n');
  insert(db, 'disruption_warning', fields, values(itime,:));
end;

fprintf(1, 'Finished successfully!\n\n');

