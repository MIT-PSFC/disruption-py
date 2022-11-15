db = set_database('d3drdb');
results = fetch(db,['select shot from disruption_warning order by shot']);
close(db)
shotlist = unique(cell2mat(results));
clearvars -except shotlist
