load('minor_disrupt_all_data_current_and_thermal_decay_unintentional_intersect_good')
shotlist=minor_disrupt_shot;
disruption_limit=cell(length(shotlist),1);
time_until_minor=cell(length(shotlist),1);
mdsconnect('atlas.gat.com');
nshots = length(shotlist);
ndigits = int2str(floor(log10(nshots)) + 1);
for i=1:length(shotlist)
    %minor_disrupt_time{i}=[];
    shot=shotlist(i);
    minor_dis_time=minor_disrupt_time{i}(1);
    fprintf(1,['Processing shot %7i  (%' ndigits 'i/%' ndigits ...
    'i  %6.2f%%)\n'], shot, i, nshots, i/nshots*100);
    
    %[shotopened, status]=mdsopen('d3d', shot);
    efittrees  = select_efit_trees(shot, '', 'DIS');
    if isempty(efittrees)
        fprintf(1, 'No disruption EFIT tree for this shot\n');
        fprintf(1, 'Using EFIT01\n');
        efittrees = {'EFIT01'};
    end
    tree=char(efittrees(end,:));
    [shotopened, status] = mdsopen(tree, shot);

    
    [efittime, status] = mdsvalue('\efit_a_eqdsk:atime');
    efittime = efittime/1000; % efit time in seconds
    
    q95 = mdsvalue('\efit_a_eqdsk:q95');
    kappa = mdsvalue('\efit_a_eqdsk:kappa');
    
    chisq = mdsvalue('\efit_a_eqdsk:chisq');
    invalid_indices = find(chisq > 50);
    kappa(invalid_indices) = NaN;
    q95(invalid_indices) =NaN;
    mdsclose;
    [~, Te_width_normalized] = get_TS_data(shot, efittime);
    safety_factor=(Te_width_normalized.^0.5).*(q95.^0.5).*(kappa.^(-1));
    safety_factor=interp1(efittime, safety_factor, efittime, 'linear');
    safety_factor=smoothdata(safety_factor,'movmean',3);
    time_need=minor_dis_time-0.1:0.01:minor_dis_time;
    dis_limit=interp1(efittime, safety_factor, time_need', 'linear');
    t_until=0.1:-0.01:0;
    indx=~isnan(dis_limit);
    disruption_limit{i}=dis_limit(indx);
    time_until_minor{i}=t_until(indx);
end