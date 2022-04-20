function [amplitude, phase, time] = ...
    get_saddle_fourier_components(shot, do_not_subtract_RMP);

if ~exist('do_not_subtract_RMP', 'var');
  do_not_subtract_RMP = 0;
elseif isempty(do_not_subtract_RMP);
  do_not_subtract_RMP = 0;
elseif do_not_subtract_RMP ~= 1;
  do_not_subtract_RMP = 0;
end;

% This routine reads in the data from the toroidal array of 8 saddle coils
% on EAST, and then does a fast Fourier transform (FFT) to determine the
% amplitude and phase of the n = 0, 1, 2, and 3 components as a function
% of time.
%
%  2016/12/05 -- Gu Xiang: write Matlab code to read in saddle coil signals
%                  and calculate Fourier component amplitudes.
%  2016/12/07 -- R. Granetz: modify Gu Xiang's code to convert it to a
%                  function, add the phase information, improve
%                  computational efficiency, and add comments
%

% -------------------------------------------------------------------
%   Check the lock mode
%   05/12/2016 by Xiang Gu
% -------------------------------------------------------------------

% Read in the saddle sensor data and the rmp currents from the MDSplus
% tree.  All the outputs have time as their 1st dimension,
% i.e. rmp(time,rmpcoil#) and saddle(time,saddlecoil#), and the time arrays
% are column vectors.

% The ordering of the rmp coils is: irmpu1 to irmpu8, then irmpl1 to irmpl8.
% The ordering of the saddle coils is: sad_pa, sad_bc, .. sad_no.

[rmptime, rmp, saddletime, saddle] = get_rmp_and_saddle_signals(shot);





% baseline_indices = find( saddletime < -5.5);
% 
% for i = 1:8;
%     acoeffs = polyfit(saddletime(baseline_indices), saddle(baseline_indices, i), 1);
%     baseline = polyval(acoeffs, saddletime); 
%     saddle(:,i) = saddle(:,i) - baseline;
% 
% end;


  
  
  
if length(saddletime) <= 2;
  amplitude = [];
  phase = [];
  time = [];
  return;
end;

% Zero_point = find( saddletime > 1.99999 & saddletime < 2.00001);
% saddle_zero_point = saddle(Zero_point, :);
% 
% for i = Zero_point:length(saddle);
%     saddle(i,:) = saddle(i,:) - saddle_zero_point;
% end;

if (do_not_subtract_RMP ~= 1);

% The saddle signals can include direct pickup from t`he RMP coils, when the
% RMP coils are active.  This direct pickup can be subtracted from the
% saddle signals, leaving just the plasma response and baseline drifts.
% The signal on each saddle sensor that is due to each RMP coil was derived
% from two RMP calibration shots, 61846 and 62053, using the Matlab routine
% "calculate_rmp_saddle_coeff_matrix.m"

  if shot > 75808 & shot < 81702;
    load('rmp_saddle_coeff_matrix_2018.mat');
    coeff_matrix(:,1) = 0;
    coeff_matrix(:,15) = 0;
  else;
    load('rmp_saddle_coeff_matrix_2017.mat');
  end;
%   load('G2.mat');
%   load('rmp_saddle_coeff_matrix_2018.mat');
%   coeff_matrix = calculate_rmp_saddle_coeff_matrix;
%   coeff_matrix(:,1) = 0;
%   coeff_matrix(:,15) = 0;
   
   

% Calculate the direct RMP contribution to the saddle sensors, and subtract
% it from the saddle signals.

%   rmp_pickup = transpose(G2 * transpose(rmp));
  rmp_pickup = transpose(coeff_matrix * transpose(rmp));
  rmp_pickup = interp1(rmptime, rmp_pickup, saddletime);
% when shot > 69635, the unit of RMP from V to A. So there is a change in
% coefficient(812).
  if shot > 69635;
%        saddle = saddle - rmp_pickup./1000;
     saddle = saddle - rmp_pickup ; 
% When shot < 56305, the coefficient of SAD is 0.02, but the coefficient of
% the calibrated shot is -0.02
  elseif shot < 56305;    
    saddle = saddle + rmp_pickup;
  else;
    saddle = saddle - rmp_pickup;
  end;
  
end;

% Zero_point = find( saddletime > 1.99999 & saddletime < 2.00001);
% saddle_zero_point = saddle(Zero_point, :);
% 
% for i = Zero_point:length(saddle);
%     saddle(i,:) = saddle(i,:) - saddle_zero_point;
% end;


%   if shot > 69622;
%       saddle = saddle - rmp_pickup./1000;
%   else;
%       saddle = saddle - rmp_pickup;







% Gushuai: Calculate fast Fourier transforms and get mode amplitudes and phases
% for i=1:size(saddle,1)
%     saddle(i,:)=saddle(i,:)-mean(saddle(i,:),2);
% end

output = fft(saddle, [], 2);  % Take FFT along 2nd dimension (phi)

amplitude = abs(output)/size(output,2);
amplitude(:, 2:end) = amplitude(:, 2:end) * 2; 

phase = atan2(imag(output), real(output));

% Only want n = 0, 1, 2, 3 (first four Fourier modes)

amplitude = amplitude(:,1:4);
phase = phase(:, 1:4);
time = saddletime;

% % the first feature 

% indices = find(time>=0.0032);
% time_first_feature = time(indices);
% for i=time_first_feature;
%     indices_interval = find(time>=i-0.0032 & time<i);
    
    


% dt = mean(diff(time));
% npts = round(0.0032/dt);
% M = movmean(amplitude, [npts, 0], 1) 
% f1 = amplitude-M;
% 
% for i = 1:length(time);
%     time(i)
% move_time = find(time)


return;
end


% figure;plot(time, abs(gradient(smooth(amplitude(:,2),501))));xlim([2,4.5]);
