function n_mode_fit_fun = n1_mode_fit_function(x, xdata)

% This function is used to fit a cosine function to the measured poloidal
% magnetic field fluctuations at different toroidal locations for a SINGLE
% time t.
%
% NOTE: This code is the same thing as n_mode_fit_function.m in
% /home/tinguely/Disruptions/Code.
%
% Author: Alex Tinguely 2015-09-11
%
% Inputs:
%   x = array of parameters to be fit
%   x(1) = amplitude of the fluctuations of the B-field in time (T/s)
%   x(2) = phase of the B-field fluctuations (rad)

n = 1; % here we are considering the n=1 MHD mode, but it could be changed

n_mode_fit_fun = x(1)*cos(n*xdata - x(2)); 
% d/dt (B_poloidal_fluctuation) = amplitude*cos(n*angles - phase)

end