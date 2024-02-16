db = set_database('logbook');
result = fetch(db, ['select dbkey, p_rad from disruption_warning ' ...
  'order by shot, time']);
dbkey = int32(cell2mat(result(:,1)));
p_rad = cell2mat(result(:,2));

p_rad = p_rad * 4.5; % Factor of 4.5 comes from cross-calibration with
                     % 2pi_foil during flattop times of non-disruptive
                     % shots, excluding times for
                     % which p_rad (uncalibrated) <= 1.e5 W

p_rad_calibrated = num2cell(p_rad);

for ikey = 1:length(dbkey);
  update(db, 'disruption_warning', {'p_rad'}, p_rad_calibrated(ikey), ...
    ['where dbkey = ' num2str(dbkey(ikey),'%i')]);
end;
close(db);
