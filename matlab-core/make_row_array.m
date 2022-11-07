function [row_array, original_orientation] = make_row_array(array)

% This function turns all arrays of lenght n (1xn row arrays or nx1 column
% arrays) into 1xn row arrays. This is used to make calculations easier.
%
% Input:
%   array = 1xn or nx1 array
%
% Outputs:
%   row_array = 1xn array
%   original_orientation = 0 if originally column array, 1 if originally row array
%
% Author: Alex Tinguely, 2015-09-30

dims = size(array); % dimensions of array, [1,n] for 1xn row array OR [n,1] for nx1 column array

if dims(1) <= dims(2) % i.e. 1 < n or row array

    row_array = array; % keep it as a row array
    original_orientation = 1; % "true", originally was a row array
    
else % if dims(1) > dims(2), i.e. n > 1 or column array
    
    row_array = transpose(array); % transpose to row array
    original_orientation = 0; % "false", originally was NOT a row array
    
end

end