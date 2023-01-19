function [data,time,z] = get_TS_cmod(shot)

[~, status] = mdsopen('electrons',shot);

% Read in Thomson core temperature data, which is a 2-D array, with the
% dependent dimensions being time and z (vertical coordinate)

node = '.yag_new.results.profiles';

[~,statA] = mdsvalue([node ':te_rz']);
[~,statB] = mdsvalue(['dim_of(' node ':te_rz)']);
[~,statC] = mdsvalue([node ':z_sorted']);

if (mod(statA,2) == 1 & mod(statB,2) == 1 & mod(statC,2) == 1)
    data = mdsvalue([node ':te_rz']);
    time = mdsvalue(['dim_of(' node ':te_rz)']);
    z = mdsvalue([node ':z_sorted']);
    mdsclose;
else;
    data = NaN;
    time = NaN;
    z = NaN;
    mdsclose;
    return;
end;
