db = set_database('logbook');
result = fetch(db, ['select dbkey, radiated_fraction from ' ...
  'disruption_warning order by shot, time']);
dbkey = int32(cell2mat(result(:,1)));
radiated_fraction = cell2mat(result(:,2));

radiated_fraction = radiated_fraction * 4.5; % Factor of 4.5 comes from
                     % cross-calibration with
                     % 2pi_foil during flattop times of non-disruptive
                     % shots, excluding times for
                     % which p_rad (uncalibrated) <= 1.e5 W

radiated_fraction_calibrated = num2cell(radiated_fraction);

nrecords = length(dbkey);
for ikey = 1:nrecords;
  fprintf(1,['Processing record # %6i/%6i (%8.4f%%)\n'], ikey, ...
  nrecords, ikey/nrecords*100);

  update(db, 'disruption_warning', {'radiated_fraction'}, ...
    radiated_fraction_calibrated(ikey), ...
    ['where dbkey = ' num2str(dbkey(ikey),'%i')]);
end;
close(db);
