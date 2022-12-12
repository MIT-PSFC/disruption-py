function [ vals, times, ier_post ] = d_geom_test2(sig_name, data)

        diff_psi_frac_thresh_is_bdy = 1e-4;
        diff_psi_frac_thresh_is_ext_sep = 0.1;

        diff_psi_frac_thresh_is_almost_perfect_dnd = 0.01;
        % diff_psi_frac_thresh_is_perfect_dnd = diff_psi_frac_thresh_is_bdy;


        ier_post = 0;

        if ~isfield(data,'efsg1psi') || isnan(data.ier_efsg1psi) || data.ier_efsg1psi>0
            vals=nan;
            times=nan;
            ier_post=1;
            return
        end
        if ~isfield(data,'efsg2psi') || isnan(data.ier_efsg2psi) || data.ier_efsg2psi>0
            vals=nan;
            times=nan;
            ier_post=1;
            return
        end
        if ~isfield(data,'efsg1zx') || isnan(data.ier_efsg1zx) || data.ier_efsg1zx>0
            vals=nan;
            times=nan;
            ier_post=1;
            return
        end
%         if ~isfield(data,'efsg2zx') || isnan(data.ier_efsg2zx) || data.ier_efsg2zx>0
%             vals=nan;
%             times=nan;
%             ier_post=1;
%             return
%         end
        if ~isfield(data,'efspsibdy1') || isnan(data.ier_efspsibdy1) || data.ier_efspsibdy1>0
            vals=nan;
            times=nan;
            ier_post=1;
            return
        end
        if ~isfield(data,'efspsimag') || isnan(data.ier_efspsibdy1) || data.ier_efspsibdy1>0
            vals=nan;
            times=nan;
            ier_post=1;
            return
        end


        g2_off = sum(abs(data.efsg2psi)) < diff_psi_frac_thresh_is_bdy;

        diff_norm_denom = data.efspsibdy1 - interp1(data.efspsimag_t,data.efspsimag,data.efspsibdy1_t);

        diff_g1 = data.efsg1psi - data.efspsibdy1;
        diff_g1 = diff_g1 ./ diff_norm_denom;
        diff_g1 = abs(diff_g1);

        g1_is_bdy = diff_g1 < diff_psi_frac_thresh_is_bdy;

        if ~g2_off
            diff_g2 = data.efsg2psi - data.efspsibdy1;
            diff_g2 = (data.efsg2psi - data.efspsibdy1) ./ diff_norm_denom;
            g2_is_bdy = abs(diff_g2) < diff_psi_frac_thresh_is_bdy;
        else
            diff_g2 = nan;
            g2_is_bdy = 0;
        end


        if g2_off       % Need to figure out if USN or LSN:
            % The g1 definitions is dependent on which is the current
            % limiting surface. CHeck g1's x-pt's Z to check if its the UN
            % or the LN:
            g1_is_ln = data.efsg1zx < -diff_psi_frac_thresh_is_bdy; % If the null is below zero, then LN
        else
            % If g1 and g2 both
            g1_is_ln = 1; %
        end

        is_lim = ~g1_is_bdy & ~g2_is_bdy;

        is_dnd_up_or_usn = ~is_lim & ((g1_is_bdy & ~g1_is_ln) | (g2_is_bdy & g1_is_ln));
        is_dnd_down_or_lsn = ~is_lim & ((g1_is_bdy & g1_is_ln) | (g2_is_bdy & ~g1_is_ln));

        if g2_off
            is_dnd=0;
            is_dnd_perfect = 0;
            is_dnd_almost_perfect=0;
        else
            is_dnd = ~is_lim & (diff_g1 < diff_psi_frac_thresh_is_ext_sep) & (diff_g2 < diff_psi_frac_thresh_is_ext_sep);
            is_dnd_perfect = g1_is_bdy & g2_is_bdy;
            is_dnd_almost_perfect = is_dnd_perfect & (diff_g1 < diff_psi_frac_thresh_is_almost_perfect_dnd) & (diff_g2 < diff_psi_frac_thresh_is_almost_perfect_dnd);
        end

        is_sn = ~is_lim & ~is_dnd;

        is_lsn = is_sn & is_dnd_down_or_lsn;
        is_usn = is_sn & is_dnd_up_or_usn;

        is_dnd_down = ~is_lsn & is_dnd_down_or_lsn;
        is_dnd_up = ~is_usn & is_dnd_up_or_usn;

        almost_exact_dnd_slight_down = is_dnd_down & is_dnd_almost_perfect;
        almost_exact_dnd_slight_up = is_dnd_up & is_dnd_almost_perfect;

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
        geom = ~is_lim;
        geom = geom(:) .* ( 1*is_lsn(:) + 2*is_dnd_down(:) + 2.4*almost_exact_dnd_slight_down(:) + 2.5*is_dnd_perfect(:) + 2.6*almost_exact_dnd_slight_up(:) + 3*is_dnd_up(:) + 4*is_usn(:) );

        vals = geom;
        times = data.efsg1psi_t;


        return
end