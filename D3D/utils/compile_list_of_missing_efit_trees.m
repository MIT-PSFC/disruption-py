load('~/disruption_warning_database/d3d_shot_duration_disrupt.mat');

list_of_missing_efit_trees = [];

db = set_database('code_rundb');

for ishot = 1:length(shotlist);
  shot = shotlist(ishot);
  efit_tree_list = fetch(db, ['select tree from plasmas where ' ...
    'shot = ' num2str(shot) ' and ' ...
    'run_by = ''granetzr'' and ' ...
    'runtag = ''DIS'' and ' ...
    'run_id = ' num2str(shot) ' and ' ...
    'deleted = 0 order by idx']);
  if isempty(efit_tree_list);
    list_of_missing_efit_trees = [list_of_missing_efit_trees, shot];
  end;
end;

close(db)
save('list_of_missing_efit_trees.mat', 'list_of_missing_efit_trees');
