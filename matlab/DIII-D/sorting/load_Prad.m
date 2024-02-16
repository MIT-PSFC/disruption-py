function data = load_Prad(pulse,fan,smoothing_window,EFIT)
%LOAD_PRAD Retrieves DIII-D radiation data from the bolometer MDS+ tree
%   Author: Kevin Montes
%   Date: March 2019
%	Inputs:
%		- pulse: shot number for which to pull data
%		- fan: configuration of channels to use ('lower','upper',or 'custom')
%		- smoothing_window: time window [in ms] to use for raw signal filter

% Add path for the 
addpath('/fusion/projects/disruption_warning/peaking_factors_d3d/Physics-based_indicators/DIAG_parameterization');

% Focus on either top or bottom fan (two separate peaking factors)
if strcmp(fan,'upper')
  fan_chans = 1:24;
elseif strcmp(fan,'lower')
  fan_chans = 25:48;
elseif strcmp(fan,'custom')
  fan_chans = [3 4 5 6 7 8 9 12 14 15 16 22]+24; % 1st choice (heavily cover divertor and core)
end

% Get causally-filtered bolometer data using Bob's routine 'getbolo.m',
% a translation of Tony Leonard's 'getbolo.pro' IDL routine.
% This may require extra processing, but not sure until talking to Bob
a = getbolo_new(pulse, smoothing_window);
b = powers_new(a);

% Sometimes the bolo data is garbage.  Check the 'ier' subfields to
% determine this
ier = 0;
ch_avail = [];
x = NaN(length(EFIT.time),length(fan_chans));
z = [];
brightness = [];
data.power = [];
if ~strcmp(fan,'custom')
	for i = 1:length(fan_chans);
	  ichan = fan_chans(i);
	  if a.chan(ichan).ier == 0;
		ch_avail = [ch_avail,ichan];
	  end;
	  x(:,i) = a.chan(ichan).Z + tan(a.chan(ichan).angle*pi/180)*(EFIT.rmaxis-a.chan(ichan).R);
	  b.chan(ichan).chanpwr(b.chan(ichan).chanpwr<0) = 0;
	  b.chan(ichan).brightness(b.chan(ichan).brightness<0) = 0;
	  z = [z;b.chan(ichan).chanpwr];
	  brightness = [brightness; b.chan(ichan).brightness];
	end;
	data.power = z;
else
	lower_fan_chans = 25:48;
	j = 0;
	for i = 1:length(lower_fan_chans);
		data.power = [data.power;b.chan(lower_fan_chans(i)).chanpwr];
		if ismember(lower_fan_chans(i),fan_chans)
			j = j+1;
			ichan = fan_chans(j);
			if a.chan(ichan).ier == 0;
				ch_avail = [ch_avail,ichan];
			end;
			x(:,j) = a.chan(ichan).Z + tan(a.chan(ichan).angle*pi/180)*(EFIT.rmaxis-a.chan(ichan).R);
			b.chan(ichan).chanpwr(b.chan(ichan).chanpwr<0) = 0;
			b.chan(ichan).brightness(b.chan(ichan).brightness<0) = 0;
			z = [z;b.chan(ichan).chanpwr];
			brightness = [brightness; b.chan(ichan).brightness];
		end
	end;
end
data.x = x';
data.z = z;
data.t = a.rawtime;
data.xtime = EFIT.time';
data.ch_avail = ch_avail;
data.brightness = brightness;
