function Prad_peaking=get_Prad_peaking(shot, timebase)

if (size(timebase,2) > 1)
  timebase = transpose(timebase);
end

Prad_peaking = NaN(size(timebase));

[~, status] = mdsopen('cmod',shot);
if (mod(status,2)==1)
   [~,statA] = mdsvalue('\efit_aeqdsk:rmagx');
   [~,statB] = mdsvalue('\efit_aeqdsk:aminor');
   [~,statC] = mdsvalue('\efit_aeqdsk:zmagx');
   a=(mod(statA,2) == 1 & mod(statB,2) == 1 & mod(statC,2) == 1);
   if a
    r0=0.01*mdsvalue('\efit_aeqdsk:rmagx');
    z0=0.01*mdsvalue('\efit_aeqdsk:zmagx');
    aminor=mdsvalue('\efit_aeqdsk:aminor');
    efittime=mdsvalue('dim_of(\efit_aeqdsk:aminor)');
   else
       mdsclose;
       return
   end
   
else
    mdsclose;
    return
    
mdsclose;

end

if length(efittime)<2
   return
end

[~, status] = mdsopen('SPECTROSCOPY',shot);
if (mod(status,2)==0)
    mdsclose;
    return
end

[t_axa,status1] = mdsvalue('dim_of(\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXA:BRIGHT, 1)');
[t_axj,status2] = mdsvalue('dim_of(\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXA:BRIGHT, 1)');
[~,status3] = mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXA:BRIGHT');
[~,status4] = mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXJ:BRIGHT');

if (mod(status1,2) == 1 & mod(status2,2) == 1 & length(t_axa)>=2 & length(t_axj)>=2 & mod(status3,2) == 1 & mod(status4,2) == 1)
    r_axa=mdsvalue('dim_of(\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXA:BRIGHT)');
    r_axj=mdsvalue('dim_of(\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXJ:BRIGHT)');
    z_axa=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.AXA:Z_O');
    %z_axa=repmat(z_axa, size(timebase));
    z_axj=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.AXJ:Z_O');
    %z_axj=repmat(z_axj, size(timebase));
    bright_axa=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXA:BRIGHT');
    bright_axj=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXJ:BRIGHT');
    good_axa=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.AXA:GOOD');
    good_axa=good_axa>0;
    good_axj=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.AXJ:GOOD');
    good_axj=good_axj>0;
    r_axa=r_axa(good_axa);
    r_axj=r_axj(good_axj);
    bright_axa=bright_axa(good_axa,:);
    bright_axj=bright_axj(good_axj,:);
    aminor=interp1(efittime,aminor,timebase, 'linear');
    r0=interp1(efittime,r0,timebase, 'linear');
    z0=interp1(efittime,z0,timebase, 'linear');
    axa_interp=NaN(size(bright_axa,1),length(timebase));
    axj_interp=NaN(size(bright_axj,1),length(timebase));
    for i=1:size(bright_axa,1)
        interped=interp1(t_axa',bright_axa(i,:),timebase', 'linear');
        indx=interped<0;
        interped(indx)=nan;
        axa_interp(i,:)=interped;
    end
    for i=1:size(bright_axj,1)
        interped=interp1(t_axj',bright_axj(i,:),timebase', 'linear');
        indx=interped<0;
        interped(indx)=nan;
        axj_interp(i,:)=interped;
    end
    for i=1:length(timebase)
        axa_dist=sqrt((r_axa-r0(i)).^2+(z0(i)-z_axa)^2);
        axj_dist=sqrt((r_axj-r0(i)).^2+(z0(i)-z_axj)^2);
        axa_core_index=axa_dist<0.2*aminor(i);
        axj_core_index=axj_dist<0.2*aminor(i);
        core_radiation=[axa_interp(axa_core_index,i); axj_interp(axj_core_index,i)];
        all_radiation=[axa_interp(:,i); axj_interp(:,i)];
        Prad_peaking(i)=mean(core_radiation,'omitnan')/mean(all_radiation,'omitnan');
    end
    
elseif (mod(status1,2) == 1 & length(t_axa)>=2  & mod(status3,2) == 1)
        r_axa=mdsvalue('dim_of(\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXA:BRIGHT)');
        z_axa=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.AXA:Z_O');
        bright_axa=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXA:BRIGHT');
        good_axa=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.AXA:GOOD');
        good_axa=good_axa>0;
        r_axa=r_axa(good_axa);
        bright_axa=bright_axa(good_axa,:);
        aminor=interp1(efittime,aminor,timebase, 'linear');
        r0=interp1(efittime,r0,timebase, 'linear');
        z0=interp1(efittime,z0,timebase, 'linear');
        axa_interp=NaN(size(bright_axa,1),length(timebase));
        for i=1:size(bright_axa,1)
            interped=interp1(t_axa',bright_axa(i,:),timebase', 'linear');
            indx=interped<0;
            interped(indx)=nan;
            axa_interp(i,:)=interped;
        end
        for i=1:length(timebase)
            axa_dist=sqrt((r_axa-r0(i)).^2+(z0(i)-z_axa)^2);
            axa_core_index=axa_dist<0.2*aminor(i);
            core_radiation=axa_interp(axa_core_index,i);
            all_radiation=axa_interp(:,i);
            Prad_peaking(i)=mean(core_radiation,'omitnan')/mean(all_radiation,'omitnan');
        end
elseif (mod(status2,2) == 1 & length(t_axj)>=2 & mod(status4,2) == 1)
    r_axj=mdsvalue('dim_of(\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXJ:BRIGHT)');
    z_axj=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.AXJ:Z_O');
    bright_axj=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXJ:BRIGHT');
    good_axj=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.AXJ:GOOD');
    good_axj=good_axj>0;
    r_axj=r_axj(good_axj);
    bright_axj=bright_axj(good_axj,:);
    aminor=interp1(efittime,aminor,timebase, 'linear');
    r0=interp1(efittime,r0,timebase, 'linear');
    z0=interp1(efittime,z0,timebase, 'linear');
    axj_interp=NaN(size(bright_axj,1),length(timebase));
    for i=1:size(bright_axj,1)
        interped=interp1(t_axj',bright_axj(i,:),timebase', 'linear');
        indx=interped<0;
        interped(indx)=nan;
        axj_interp(i,:)=interped;
    end
    for i=1:length(timebase)
        axj_dist=sqrt((r_axj-r0(i)).^2+(z0(i)-z_axj)^2);
        axj_core_index=axj_dist<0.2*aminor(i);
        core_radiation=axj_interp(axj_core_index,i);
        all_radiation=axj_interp(:,i);
        Prad_peaking(i)=mean(core_radiation,'omitnan')/mean(all_radiation,'omitnan');
    end
else
    mdsclose;
    return
end

end
    
    
    
    
    
