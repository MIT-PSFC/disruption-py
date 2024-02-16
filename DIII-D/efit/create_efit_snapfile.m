function create_efit_snapfile(shot);

% There are two EFIT input namelist parameters, "iavem" and "iavev", that
% specify the averaging times for the magnetics data and loop voltage data
% respectively.  These are usually set to 5 (i.e. window width of +/- 5 ms)
% and 10 (+/- 10 ms) respectively.  For our disruption warning purposes, we
% want to use the actual EFIT snap file that was used for the default EFIT
% (EFIT01), except that want to change both the "iavem" and "iavev"
% parameters to +/- 2 ms.  We do this by using Matlab's string parsing
% routines.  (I initially set iavem and iavev to +/- 1 ms, but I found that
% EFIT failed on about 250 shots.  I then found that increasing the
% smoothing windows to +/- 2 ms reduced the EFIT failure rate.)

if ~exist('shot');
  shot = str2num(getenv('SHOT'));
end;

% The text node 'namelist' (at the top of every EFIT tree), contains an
% exact copy of the EFIT snap file used to produce the data in each EFIT
% tree.

mdsconnect('atlas.gat.com');
[~, status] = mdsopen('efit01', shot);

% If the EFIT01 tree exists, then read the namelist from it.  If the EFIT01
% tree does not exist, or if the namelist data does not contain the
% necessary namelist ("efitin"), then read in the namelist from the default
% EFIT snap file.  This must be done line-by-line, accumulating the lines
% into a single long string.

if (mod(status,2)==1);
  namelist = mdsvalue('namelist');
  mdsclose;
  if isempty(strfind(lower(namelist), '&efitin'));
    status = 0;
  end;
end;

if (mod(status,2)==0);
  fileid = fopen( ...
    '/fusion/projects/codes/efit/support_files/snapfiles/efit_snap.dat_jt')
  namelist='';
  end_of_file = 0;
  while ~end_of_file;
    textline = fgets(fileid);
    if (textline == -1);
      end_of_file = 1;
    else;
      namelist = [namelist textline];
    end;
  end;
  fclose(fileid);
end;

% Change the value for "iavem" to +/- 2 ms

indx = strfind(lower(namelist), 'iavem');
jndx = strfind(lower(namelist(indx(1):end)), '=');
start_of_field_indx = indx(1) + jndx(1);
letter_indices = find(isstrprop(namelist(start_of_field_indx:end), 'alpha'));
cntrl_indices  = find(isstrprop(namelist(start_of_field_indx:end), 'cntrl'));
jndx = min([letter_indices, cntrl_indices]);
end_of_field_indx = start_of_field_indx + jndx - 2;
namelist(start_of_field_indx:start_of_field_indx) = '2';
for i = start_of_field_indx+1:end_of_field_indx;
 namelist(i:i) = ' ';
end;

% Change the value for "iavev" to +/- 2 ms

indx = strfind(lower(namelist), 'iavev');
jndx = strfind(lower(namelist(indx(1):end)), '=');
start_of_field_indx = indx(1) + jndx(1);
letter_indices = find(isstrprop(namelist(start_of_field_indx:end), 'alpha'));
cntrl_indices  = find(isstrprop(namelist(start_of_field_indx:end), 'cntrl'));
jndx = min([letter_indices, cntrl_indices]);
end_of_field_indx = start_of_field_indx + jndx - 2;
namelist(start_of_field_indx:start_of_field_indx) = '2';
for i = start_of_field_indx+1:end_of_field_indx;
 namelist(i:i) = ' ';
end;

% Change the comment (first line) to describe our changes to the EFIT snap
% input file.

cntrl_indices  = find(isstrprop(namelist, 'cntrl'));  % find new line char
end_of_line_indx = cntrl_indices(1);
namelist = [ ...
 ' EFIT_SNAP.DAT_DISRUPT  (iavem and iavev set to +/- 2 ms by R. Granetz)', ...
 namelist(end_of_line_indx:end)];

% Add a blank line at end to exactly match actual efit snap files.

namelist = [namelist char(10)];

% Now write the namelist string out to a new EFIT snap file in the
% appropriate directory.

dirpath = ['shot' num2str(shot,'%i') '/'];
fileid = fopen([dirpath 'efit_snap.dat_disrupt'], 'w');
fprintf(fileid, '%s', namelist);
fclose(fileid);
