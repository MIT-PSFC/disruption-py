function data = load_bolometer(pulse,fan,varargin)
%LOAD_THOMSON Retrieves DIII-D radiation data from the bolometer MDS+ tree
%   Author: Kevin Montes
%   Date: March 2019

measure = 'chanpwr'; % set to either 'brightness' or 'chanpwr'

% Focus on either top or bottom fan (two separate peaking factors)
if strcmp(fan,'upper')
  fan_chans = 1:24;
elseif strcmp(fan,'lower')
  fan_chans = 25:48;
elseif strcmp(fan,'custom')
  %fan_chans = [2:2:24]+24;
  fan_chans = [3 4 5 6 7 8 9 12 14 15 16 22]+24; % 1st choice (heavily cover divertor and core)
  %fan_chans = [2 4 6 8 10 12 14 16 18 19 22 24]+24; % 2nd choice (sparse -> noisier and different baselines)
end

% Get causally-filtered bolometer data using Bob's routine 'getbolo.m',
% a translation of Tony Leonard's 'getbolo.pro' IDL routine.
% This may require extra processing, but not sure until talking to Bob

smoothing_window = 0.040; % use 40 ms causal smoothing window
a = getbolo_new(pulse, smoothing_window*1.e3);
b = powers_new(a);
% Sometimes the bolo data is garbage.  Check the 'ier' subfields to
% determine this
ier = 0;
ch_avail = [];
% grab EFIT data for this particular shot
if length(varargin) > 0
    trees_to_try = varargin{1};
    EFIT = load_efit(pulse,trees_to_try);
else
    EFIT = load_efit(pulse); 
end
x = NaN(length(EFIT.time),length(fan_chans));
z = [];
for i = 1:length(fan_chans);
  ichan = fan_chans(i);
  if a.chan(ichan).ier == 0;
    ch_avail = [ch_avail,ichan];
  end;
  x(:,i) = a.chan(ichan).Z + tan(a.chan(ichan).angle*pi/180)*(EFIT.rmaxis-a.chan(ichan).R);
  chan_data = b.chan(ichan).(measure);
  chan_data(chan_data<0) = 0;
  z = [z;chan_data];
end;
data.x = x';
data.z = z;
data.t = a.rawtime;
data.xtime = EFIT.time';
data.ch_avail = ch_avail;

%{
% Assign channel configurations for different plasma shapes
% Channels 1-24 (U01-U24) and 25-48 (L01-L24)
configs = struct('SNB',struct,'SNT',struct,'DN',struct);
configs.SNB.core = [38:43]; configs.SNB.edge = [33:35];
configs.SNT.core = [9:14]; configs.SNT.edge = [17:20];
configs.DN.core = []; configs.DN.edge = [];

% Loop through each time and populate the peaking factors
indices = find(bolo.rawtime>min(EFIT.time) & bolo.rawtime<max(EFIT.time));
peaking_factor = NaN(size(indices));
for i = 1:length(indices)
    t = bolo.rawtime(indices(i));
    [~,efit_indx] = min(abs(EFIT.time-t));
    c = find(strcmp(fieldnames(configs),EFIT.limloc(efit_indx)));
    if ~isempty(c)
        % Plasma is single null or double null choose one fan or other, 
        % depending on where the X-point is
        channels = configs.(string(EFIT.limloc(efit_indx)));
    else
        % Not sure what to do here? For now, assume plasma is double null
        channels = configs.('DN');
    end
    % Now, loop through channels and take average value of core/edge
    core = []; edge = [];
    for j=1:length(channels.core)
        core = [core,bolo.chan(channels.core(j)).pwr(indices(i))];
    end
    for j=1:length(channels.edge)
        edge = [edge,bolo.chan(channels.edge(j)).pwr(indices(i))];
    end
    peaking_factor(i) = sum(core,'omitnan')/sum(edge,'omitnan');
end
%}
