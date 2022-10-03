db = set_database('logbook');
exec(db, 'alter table disruption_warning add dipprog_dt float');
close(db);
