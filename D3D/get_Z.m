function [Z_error, Z_error_normalized Z_prog, Z_cur, Z_cur_normalized] = ...
    get_Z_parameters(shot, timebase);

% This script was adapted from get_Z_parameters.m on C-Mod.
%
% The purpose of this script is to read in the values of Z_error and Z_prog
% from the plasma control system (PCS). Z_prog is the programmed vertical
% position of the plasma current centroid, and Z_error is the difference
% between the actual position and that requested (Z_error = Z_cur -
% Z_prog). Thus, the actual (estimated) position, Z_cur, can be calculated.
% We do not yet know the programmed Z_prog.
% These values are then linearly interpolated over the given timebase.
%
% Inputs:
%   shot = shot number
%   timebase = array of desired time values [s]
%
% Outputs:
%%   Z_prog = programmed (i.e. requested) vertical position of the current
%     centroid [m]
%   Z_cur = actual vertical position of the current centroid, calculated
%         = Z_cur - Z_prog [m]
%
% Original Author: Alex Tinguely 2015-09-09
% Rewritten by Robert Granetz 2016/01/29
% Updated by Alex Tinguely 2016-05-02
% Updated by Robert Granetz 2017/12/22 to use Ip from magnetics tree instead
%   of the PCS tree for shots prior to 2015.  This is because a factor
%   required to convert the PCS Ip signal to units of amperes is apparently
%   not available in the PCS tree for shots prior to the 2015 run campaign
% Updated by Cristina Rea 2019/08/23 for DIII-D.


% The input array, "timebase", can be either a row vector (1 x n) or a
% column vector (n x 1).  We desire to have the output arrays match the
% shape of "timebase".  However, Matlab calls to the routine "mdsvalue"
% return column vectors for 1-D signals, so it is simpler to work with
% column vectors within this routine, and then, if necessary, convert the
% outputs to row vectors just before exiting this routine.  So the first
% step is to create a copy of "timebase" that is guaranteed to be a column
% vector.

if (size(timebase,1) > 1);
  timebase_column = timebase;
else;
  timebase_column = transpose(timebase);
end;

% Initialize all output arrays to NaN (Not-a-Number) column vectors

Z_cur = NaN(length(timebase),1);

