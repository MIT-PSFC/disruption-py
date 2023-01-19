% <alessandro.pau@protonmail.ch>
% Department of Electrical and Electronic Engineering of University of Cagliari
% Created: 2015
% Modified: 2018

function data = bolo_coord_conv(data,RMAG)

F1 = [4.945, 0.336]; 
F2 = [4.993, 0.336]; 
F1_MAG_R = F1(1)-RMAG.signal;
F2_MAG_R = F2(1)-RMAG.signal;
Mat = [repmat(F1_MAG_R,1,8),repmat(F2_MAG_R,1,12),repmat(F1_MAG_R,1,4)];
Mat = Mat(:,data.ch_avail);
data.x = bsxfun(@times,Mat,tan(deg2rad(180-data.x)))+F1(2);  
data.x = data.x';
data.x_old = data.x;
data.xtime = RMAG.time';

function rad = deg2rad(deg)
rad = (pi/180).*deg;