function Profiles = get_peaking_factors_d3d(pulse,t_window,bolo_fan,varargin)

%-----------------------------------------------------------------------
% This function calculates peaking factors for the shot/pulse numbers
% given by the user during the times that are specified in t_window.
% It then plots the peaking factors along with the raw data and a few 
% other relevant signals in the PF_scope routine. A Matlab structure 
% with all the data and peaking factor information is returned.
%
% Inputs:
%	- pulse: shot number, or an Nx1 array of shot numbers for N shots
%	- t_window: an Nx2 array with lower and upper time limits for N shots
%	- bolo_fan: a string specifying the bolometer fan to use for the
%		    P_rad peaking factors (can be either 'lower' or 'upper')
% 	- efit_tree (optional): can specify a string w/ name of EFIT
%				tree to use (for example, 'EFIT03')
% Outputs:
%	- Profiles: an N by 4 structure containing Te, ne, P_rad (core 
% 		    vs. all), and P_rad (divertor vs all) information
% 		    for the N shots specifed in 'pulse'
%
% Author: Kevin Montes
% Date: June 24, 2019
%-----------------------------------------------------------------------

% Add libraries & paths
% ---------------------
addpath('/fusion/projects/disruption_warning/peaking_factors_d3d/Physics-based_indicators/DIAG_parameterization');
addpath('/fusion/projects/disruption_warning/peaking_factors_d3d/Physics-based_indicators');

% If needed, change data type of pulse to be consistent with t_window
if ~isa(pulse,'double')
    pulse = cast(pulse,'double');
end

% Diagnostics array definition -> want Te & ne from Thomson scattering as well as P_rad from foil bolometer
Diag_array = {'TS_{TE}';'TS_{NE}';'BOLO_{CVA}';'BOLO_{X-DIV}'};

% Check pulses and t_window size
if numel(pulse) ~= size(t_window,1);
    disp('The number of pulses has to be consistent with t_window!')
    return
end
fprintf('\n  Pulse %d: profile calculation started \n',pulse)

% Call routine for profile generation (remember to set properly save_mode - the default setting is reading the user)    
pulse_catch = [];
Profiles = cell(numel(pulse),numel(Diag_array));
for jj = 1 :numel(pulse) % can change to parfor loop to parallelize if doing long list of shots
    timebase = [t_window(1): 0.005: t_window(end)]; % sampling period of 5 ms
    try
	if length(varargin) > 0
	    efit_tree = {varargin{1}};
            Profiles(jj,:) = DIS_tool_profiles_gen(pulse(jj),timebase,Diag_array,bolo_fan,efit_tree);
	else
	    Profiles(jj,:) = DIS_tool_profiles_gen(pulse(jj),timebase,Diag_array,bolo_fan);
	end
    catch
        pulse_catch = [pulse_catch; pulse(jj)]; %#ok
        disp(['Profile gen not completed for pulse ',num2str(pulse(jj))])
    end  
    
    if mod(round(jj/numel(pulse)*100),5) == 0
        disp([num2str(jj/numel(pulse)*100),'% done ...'])
    end
    %if length(varargin) > 0
	%trees_to_try = {varargin{1}};
	%PF_scope(pulse(jj),Profiles,Diag_array,pulse,'auto',trees_to_try)
    %else
	%PF_scope(pulse(jj),Profiles,Diag_array,pulse,'auto')
    %end

end

clear jj 

fprintf('\n profile calculation completed \n')

% Save all data
%save('PF_SNT_2016_core30_vs_all.mat')
