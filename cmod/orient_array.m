function reoriented_array = orient_array(array, original_orientation)

% This function works in conjunction with make_row_array.m. If the original
% orientation of the timebase array in make_row_array() is 0, so a column
% array, this script turns the given array into a column array. If
% original_orientation is 1, we turn the array into a row array.
%
% Inputs:
%   array = 1xn row array or nx1 column array
%   original_orientation = 0 if originally column array, 1 if originally row array
%
% Output:
%   reoriented_array = row or column array depending on value of original_orientation
%
% Author: Alex Tinguely, 2015-09-30

dims = size(array);

if (original_orientation == 0 && dims(1) <= dims(2)) || (original_orientation == 1 && dims(1) > dims(2))
% if we want a column array and array is a row array OR if we want a row array and array is a column array  
    reoriented_array = transpose(array);
    
else % if (original_orientation == 0 && dims(1) > dims(2)) || (original_orientation == 1 && dims(1) <= dims(2))
% want a column array and array is a column array OR want a row array and array is a row array
    reoriented_array = array;
    
end

end