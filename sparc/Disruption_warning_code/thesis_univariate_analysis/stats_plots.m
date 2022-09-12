% Individual univariate histograms for each machine
toks = {'cmod','d3d','east'};
t_class = {0.04,0.35,0.1};
var_to_hist = 'beta_frac';

% Scatter plot/histogram with two vars
scatter_hist = {{'d3d'},{'Greenwald_fraction','beta_frac'}};

%------------------------------------------------------------------------
tbl_all = table('Size',[0,3],'VariableNames',{var_to_hist,'Class','Tokamak'}, ...
			'VariableTypes',{'double','string','string'});
for i=2:2%length(toks)
	load([toks{i} '_sql_db.mat']);
	if strcmp(toks{i},'east') & strcmp(var_to_hist,'radiated_fraction')
		var_to_plot = 'rad_input_frac';
	elseif strcmp(var_to_hist,'ip_error_frac')
		ip_prog = vars.ip - vars.ip_error;
		vars.ip_error_frac = vars.ip_error./ip_prog;
	elseif strcmp(var_to_hist,'beta_frac')
		vars.beta_frac = vars.beta_n./(4*vars.li);
	end
	% Find set of intentional disruptions and hardware failure shots to exclude
	to_exclude = intersect(vars.intentional_disruption,vars.power_supply_railed);
	if isfield(vars,'other_hardware_failure')
		to_exclude = intersect(to_exclude,vars.other_hardware_failure);
	end
	% Divide indices into disruptive and non-disruptive
	near_disruption = find(time_until_disrupt<t_class{i});
	flattop_near_disrupt = intersect(indices_flattop_disrupt_in_flattop,near_disruption);
	flattop_far_disrupt = setdiff(indices_flattop_disrupt_in_flattop,near_disruption);
	far_from_d_indx = setdiff(flattop_far_disrupt,to_exclude);
	d_indx = setdiff(flattop_near_disrupt,to_exclude);
	nd_indx = setdiff(indices_flattop_no_disrupt,to_exclude);

	% Plot histogram
	figure()
	p = prctile(vars.(var_to_hist),[2,96]);
	edges = linspace(p(1),p(2));
	histogram(vars.(var_to_hist)(nd_indx),edges,'normalization', 'probability', ... 
		'facecolor','blue','facealpha',0.3)
	title([toks{i}])
	xlabel(var_to_hist)
	ylabel('Probability')
	hold on;
	histogram(vars.(var_to_hist)(far_from_d_indx),edges,'normalization', 'probability', ... 
		'facecolor','green','facealpha',0.3)
	histogram(vars.(var_to_hist)(d_indx),edges,'normalization', 'probability', ... 
		'facecolor','red','facealpha',0.3)
	t = num2str(int32(t_class{i}*1e3));
	legend('non-disruptive', ['time until disrupt > ' t ' ms'], ...
		['time until disrupt < ' t ' ms'], 'location', 'best');

	% If selected, plot 2D scatter plot/hist
	if any(cellfun(@(x) strcmp(x,toks{i}),scatter_hist{1}))
		tbl = table('Size',[length(shot), 3],'VariableTypes', ...
			{'double','double','string'},'VariableNames', ...
			{scatter_hist{2}{1},scatter_hist{2}{2},'Class'});
		tbl.(scatter_hist{2}{1}) = vars.(scatter_hist{2}{1});
		tbl.(scatter_hist{2}{2}) = vars.(scatter_hist{2}{2});
		tbl.Class(d_indx) = 'near disrupt';
		tbl.Class(nd_indx) = 'no disrupt';
		tbl = tbl(union(d_indx,nd_indx),:);
		fig=figure();
		C = colororder; C([1,2],:) = C([2,1],:);
		colororder(fig,C);
		s = scatterhistogram(fig,tbl,scatter_hist{2}{1}, ...
				scatter_hist{2}{2},'GroupVariable','Class', ...
				'HistogramDisplayStyle','smooth','LineStyle','-')
		s.MarkerStyle = '.';
		s.MarkerSize = 10;
		s.MarkerAlpha = 0.2;
		s.XLimits = [0,15]; s.YLimits = [0,5];
		%s.Color = flip(s.Color,1);
	end
	tok_tbl = table('Size',[length(shot), 3],'VariableTypes', ...
			{'double','string','string'},'VariableNames',{var_to_hist,'Class','Tokamak'});
	tok_tbl.Tokamak(:) = toks{i};
	tok_tbl.Class(d_indx) = 'near disrupt';
	tok_tbl.Class(nd_indx) = 'no disrupt';
	tok_tbl.(var_to_hist) = vars.(var_to_hist);
	tok_tbl = tok_tbl(union(d_indx,nd_indx),:);
	tbl_all = outerjoin(tbl_all,tok_tbl,'MergeKeys', true);
end
%------------------------------------------------------------------------

tokamaks = {'C-Mod','DIII-D','EAST'};
tbl_all.Tokamak = categorical(tbl_all.Tokamak,toks,tokamaks);
tbl_all = tbl_all(find(tbl_all.(var_to_hist)~=0),:);
figure()
b = boxchart(tbl_all.Tokamak,tbl_all.(var_to_hist),'GroupByColor',tbl_all.Class, ...
			'MarkerStyle','none')
temp_color = b(1).BoxFaceColor;
b(1).BoxFaceColor = b(2).BoxFaceColor;
b(2).BoxFaceColor = temp_color;
ylabel(var_to_hist)
legend


	
