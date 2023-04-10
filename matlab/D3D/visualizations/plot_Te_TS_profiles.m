function plot_Te_TS_profiles_and_width(shot);

% Author: Robert Granetz   Dec 2017

mdsconnect('atlas.gat.com');
[shotopened, status]=mdsopen('electrons', shot);
if (mod(status,2)==1);
  [Te, status] = mdsvalue('\top.ts.blessed.core:temp'); % Read in Te(t,z)
  if (mod(status,2) == 1);                    % If successful, continue
    Z = mdsvalue('dim_of(\top.ts.blessed.core:temp, 1)');
    time = mdsvalue('dim_of(\top.ts.blessed.core:temp, 0)')/1.e3; % ms -> s
% Get rid of the last channel (#41), which is not real
    Te = Te(:, 1:end-1);
    Z = Z(1:end-1);
%   duration = end_of_shot(shot);
%   tindx = find(time > 0 & time <= duration);
    tindx = find(time > 0);
    figure(1);
    plot(Z, transpose(Te(tindx,:)));
    ylim([0,10e3]);
    set(gca, 'fontsize', 12);
    xlabel('Z [m]', 'fontsize', 14);
    ylabel('T_e [eV]','fontsize', 14);
    title(['T_e Thomson data; Shot ' num2str(shot)], 'fontsize', 16);

    Te_width_TS = NaN(length(time), 1);
    zarray = [0:.01:.9];
    for i = 1:length(tindx);
      y = Te(tindx(i),:);
      ok_indices = find(y ~= 0);
      y = y(ok_indices);
      z = Z(ok_indices);
      if (length(ok_indices) > 2);
        p = polyfit(z, transpose(y), 2);
        Te_array = polyval(p, zarray);
        [Te_max, maxindx] = max(Te_array);
        z_max = zarray(maxindx);
        Te_HM = Te_max/2;
        [~, HM_indices] = min(abs(Te_array - Te_HM));
        HM_indx = max(HM_indices);
        z_HM = zarray(HM_indx);
        if z_HM > z_max;
          Te_width_TS(tindx(i)) = z_HM - z_max;
        end;
      end;
    end;

    figure(2);
    plot(time, Te_width_TS, 'b');
    set(gca, 'fontsize', 12);
    xlabel('Time [s]', 'fontsize', 14);
    ylabel('T_e HWHM [m]', 'fontsize', 14);
    title(['T_e Thomson HWHM; Shot ' num2str(shot)], 'fontsize', 16);


    while 1;
      timeval = input('Enter time to plot profile fit (<CR> to exit): ');
      if isempty(timeval); break; end;
      figure(3);
      [~, tindx] = min(abs(time - timeval));
      timeval = time(tindx);
      plot(Z, Te(tindx,:), 'sb');
      set(gca, 'fontsize', 12)
      xlabel('Z [m]', 'fontsize', 14);
      ylabel('T_e [eV]','fontsize', 14);
      title(['T_e profile fit; Shot ' num2str(shot) ...
        ' at t = ' num2str(timeval,'%5.3f') ' s'], 'fontsize', 16);

      y = Te(tindx,:);
      ok_indices = find(y ~= 0);
      y = y(ok_indices);
      z = Z(ok_indices);
      if (length(ok_indices) > 2);
        p = polyfit(z, transpose(y), 2);
        Te_array = polyval(p, zarray);
        [Te_max, maxindx] = max(Te_array);
        z_max = zarray(maxindx);
        Te_HM = Te_max/2;
        [~, HM_indices] = min(abs(Te_array - Te_HM));
        HM_indx = max(HM_indices);
        z_HM = zarray(HM_indx);
        if z_HM > z_max;
          Te_width_TS(tindx) = z_HM - z_max;
        end;
      end;
      hold on;
      plot(zarray, Te_array, 'r');
      hold off;
    end;

   end;
  mdsclose;
end;

end
