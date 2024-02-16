%% Get IR Camera data

mdsconnect('alcdata.psfc.mit.edu')
shot = 1160617021;

[~,status] = mdsopen('spectroscopy',shot);
if mod(status,2)==1
	result = mdsvalue('\TOP.VIDEO.WIDE2:COMMENT')
else
	disp('Failed to open tree')
end
