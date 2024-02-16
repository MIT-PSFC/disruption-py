function add_column(variable_name)
    
    db = set_database('east_disruption');
    query = sprintf('ALTER TABLE disruption_warning ADD %s float', variable_name);
    exec(db, query)
    close(db)
end
