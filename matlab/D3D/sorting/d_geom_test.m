function [ vals, times, ier_post ] = d_geom_test(sig_name, data)

        seplim_thresh_div = 0.5e-2; % 1/2 cm

        drsep_thresh_absmax_DND_nobias = 1.5e-2;    %|drsep| < 1.5cm = DND w/almost no bias

        ier_post = 0;

        if ~isfield(data,'seplim') || isnan(data.ier_seplim) || data.ier_seplim>0
            vals=nan;
            times=nan;
            ier_post=1;
            return
        end

        if ~isfield(data,'drsep') || isnan(data.ier_drsep) || data.ier_drsep>0
            vals=nan;
            times=nan;
            ier_post=1;
            return
        end

        geom = data.seplim(:);

        % Large positive values mean an error:
        bad_indices = geom > 0.2; % Means bad data
        geom = geom(:) > seplim_thresh_div; % 1 if diverted, false-positive/negative if very big vals
        geom = double(geom);
        geom(bad_indices) = nan;

        % Now, what kind of diverted?
        drsep = data.drsep;

        is_lsn = abs(drsep(:)-(-0.4)) < 1e-4;
        is_lsn = double(is_lsn);

        is_usn_or_lim = abs(drsep(:)-(0.4)) < 1e-4;
        is_usn_or_lim = double(is_usn_or_lim);


        is_dn_nbl = (drsep(:) < 0) & (is_lsn(:)<=0);
        is_dn_nbl = double(is_dn_nbl);

        is_dn_pbu = (drsep(:) > 0) & ~(is_usn_or_lim(:)>0);
        is_dn_pbu = double(is_dn_pbu);


        almost_exact_dnd = abs(drsep(:)) < drsep_thresh_absmax_DND_nobias;

        almost_exact_dnd_slight_nbl = almost_exact_dnd & is_dn_nbl;
        almost_exact_dnd_slight_pbu = almost_exact_dnd & is_dn_pbu;
        exact_dnd = abs(drsep(:)) == 0;

        is_dn_nbl = is_dn_nbl & ~almost_exact_dnd;
        is_dn_pbu = is_dn_pbu & ~almost_exact_dnd;


        % Multiply by codes to designate type of div topology:
        %   geom is currently 0 or 1, now want:
        %       0 = lim
        %       1 = LSN
        %       2 = DN-B (DN neg/down biased)
        %           2.40 = In threshhold for DND no bias, slightly negative
        %           2.50 = Exact DND no bias -> drsep = 0
        %           2.60 = In threshhold for DND no bias, slightly positive
        %       3 = DB+B (DN pos/up biased)
        %       4 = USN
        geom = geom(:) .* ( 1*is_lsn(:) + 2*is_dn_nbl(:) + 2.4*almost_exact_dnd_slight_nbl(:) + 2.5*exact_dnd(:) + 2.6*almost_exact_dnd_slight_pbu(:) + 3*is_dn_pbu(:) + 4*is_usn_or_lim(:) );

        % Note: geom is still zero if limited.

        vals = geom;
        times = data.seplim_t;

        return
end