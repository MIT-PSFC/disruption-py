fileID = fopen('EAST_disruption_list.txt');

dummy_line = fgetl(fileID); % ignore comment line at top of file
dummy_line = fgetl(fileID); % ignore comment line at top of file
dummy_line = fgetl(fileID); % ignore comment line at top of file

a = fscanf(fileID, ['%i %f %f %f'],[4,inf]);

fclose(fileID);

shots = a(1,:);
tdis  = a(2,:);
ip    = a(3,:);
didt  = a(4,:);

clear a;

whos shots tdis ip didt;
