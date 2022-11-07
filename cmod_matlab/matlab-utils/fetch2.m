function result = fetch2(db, sqlquery)

% This is an attempt to mock the fetch() function in the MATLAB database
% toolbox, mostly for use on DIII-D. Its purpose is to import database data
% into MATLAB. We expect the result to be a cell array {}, which is not
% necessarily the case for the normal fetch() function.
%
% Inputs:
%   db = database. This should be defined using set_database.m for DIII-D.
%   sqlquery = string SQL query
%               "SELECT *parameter_1,...,parameter_n* FROM *table_name*
%                WHERE *some rule* ORDER BY *some ordering*"
%              example: "select shot, time from disruption_warning where
%              shot between 1150101001 and 1150105001 order by shot"
%              
% Outputs:
%   result = cell array {} of parameters requested in sqlquery
%
% Author: Alex Tinguely 2016-02-17

% If there is an error when calling the database.
if (db == 0);
  fprintf(1,'ERROR: Database is not connected.\n');
end

% I'm going to try to use information for the disruption_database.m code
% written by Bob Granetz for DIII-D
% (in /home/granetz/Collaborations/DIII-D/disruption_database). Because we
% have the Database Toolbox at C-Mod, most of these commands don't work
% here though....

% USEFUL: Result set documentation
% https://docs.oracle.com/javase/7/docs/api/java/sql/ResultSet.html

%rs.next(); % There is a cursor that points to one row of the Result Set and
           % and it starts BEFORE the first row, so we need to step
           % through using rs.next()
           
% We want to know how many rows there are in our table, corresponding to
% the number of rows in our cell array. In order to do this,
% we want to do another query using "select count(*) from..." so we parse
% our sqlquery to get the "from..." component.

sqlquery_split = regexp(sqlquery, 'from', 'split'); 
% should be CELL array of length 2, where the second index is the part of the query after "from"                                             

% NOTE: Also, we cannot have the phrase "order by ..." in our "count"
% statement, so we parse again to remove that component. However, we do
% want to include the "where ..." component. Luckily, in standard SQL query
% structure, the normal ordering is "select", "from", "where", "order by",
% so using "order" as the delimiter should work fine.

sqlquery_split2 = regexp(sqlquery_split{2}, 'order', 'split');
% this should be a cell array of length 2, where the first index includes
% the text after "from" and before "order" (so including the "where"
% component)

ps = db.prepareStatement(['select count(*) from', sqlquery_split2{1}]);
rs = ps.executeQuery();
rs.next(); % goes to first row
n_rows = rs.getObject(1); % should get us the number of rows, not sure why we use getObject, but Bob does in his code
           
% We also want to know the number of "columns" in the Result Set table,
% which will be the number of different columns in our cell array. These
% correspond to the different parameters we are querying. We should
% be able to determine this from the SQL query. We have already split up
% our SQL query, so now we split up the first part of sqlquery_split using 
% the comma as a delimiter.

sqlquery_split3 = regexp(sqlquery_split{1}, ',', 'split');
n_cols = length(sqlquery_split3); % should be number of columns

result = cell(n_rows, n_cols); % initialize output, same orientation as fetch()

ps = db.prepareStatement(sqlquery); % now execute the actual query
rs = ps.executeQuery(); % This is a Result Set, a table of data from the query
rs.next(); % go to first row of table

% Fortunately, both the column indeces in the Result Set table and matlab
% indeces start with 1. We'll read through each row and column and fill our
% cell array. We use getObject() for now, but this could be changed to
% getFloat() later, if that works.

for i = 1:n_rows
    
   for j = 1:n_cols
       
       result{i,j} = rs.getObject(j);
       rs.next(); % go to next row
       
   end
    
end

end