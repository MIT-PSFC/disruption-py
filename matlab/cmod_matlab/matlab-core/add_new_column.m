db = set_database('logbook');
exec(db, 'alter table disruption_warning add commit_hash TEXT(1');
close(db);
