function [ne_PF, Te_PF, pressure_PF] = get_peaking_factor_cmod(shot, timebase)


path('/home/jxzhu/matlab/yags', path);
addpath(genpath('/home/tinguely/Disruptions/Code/Disruption_Database'));
addpath('/home/montes/Disruption_warning_code');
addpath('/home/granetz/matlab');

if (size(timebase,1) > 1)
  timebase_column = timebase;
else
  timebase_column = transpose(timebase);
end

% Exclude all shots in the blacklist

if (shot>1120000000 && shot<1120213000) || (shot>1140000000 && shot<1140227000) || (shot>1150000000 && shot<1150610000) || (shot>1160000000 && shot<1160303000)
    ne_PF = NaN(size(timebase));
    pressure_PF = NaN(size(timebase));
    Te_PF = NaN(size(timebase));
    return
end

% Read in efit data 
[~, status] = mdsopen('cmod',shot);
if (mod(status,2)==1)
   [~,statA] = mdsvalue('\efit_aeqdsk:zmagx');
   [~,statB] = mdsvalue('\efit_aeqdsk:aminor');
   [~,statC] = mdsvalue('\efit_aeqdsk:kappa');
   a=(mod(status,2) ==1 & mod(statA,2) == 1 & mod(statB,2) == 1 & mod(statC,2) == 1);
   if a
    z0=0.01*mdsvalue('\efit_aeqdsk:zmagx');
    aminor=mdsvalue('\efit_aeqdsk:aminor');
    kappa=mdsvalue('\efit_aeqdsk:kappa');
    efittime=mdsvalue('dim_of(\efit_aeqdsk:aminor)');
    bminor=aminor.*kappa;
   else
       ne_PF = NaN(size(timebase));
       pressure_PF = NaN(size(timebase));
       Te_PF = NaN(size(timebase));
       mdsclose;
       return
   end
mdsclose;
end

if length(efittime)<2
    ne_PF = NaN(size(timebase));
       pressure_PF = NaN(size(timebase));
       Te_PF = NaN(size(timebase));
       return
end

[~, status] = mdsopen('electrons',shot);
if (mod(status,2)==0)
    ne_PF = NaN(size(timebase));
    pressure_PF = NaN(size(timebase));
    Te_PF = NaN(size(timebase));
    mdsclose;
    return
end



node = '.yag_new.results.profiles';

[~,statA] = mdsvalue([node ':ne_rz']);
[~,statB] = mdsvalue(['dim_of(' node ':ne_rz)']);
[~,statC] = mdsvalue([node ':z_sorted']);
[~,statD] = mdsvalue([node ':te_rz']);
[~,statE] = mdsvalue('\fiber_z');
[~,statF] = mdsvalue('\ts_ne');
[~,statG] = mdsvalue('\ts_te');

nlnum=4;

if (mod(statA,2) == 1 && mod(statB,2) == 1 && mod(statC,2) == 1 && mod(statE,2) == 1 && mod(statF,2) == 1)
    try
    [nlts1, nlts2, nltci1, nltci2,~,~]=compare_ts_tci(shot, nlnum);
    catch
        ne_PF = NaN(size(timebase));
        pressure_PF = NaN(size(timebase));
        mdsopen('electrons',shot);
     if mod(statD,2) == 1&& mod(statG,2) == 1
        TS_te = mdsvalue([node ':te_rz'])*1000*11600;
        tets_edge = mdsvalue('\ts_te')*11600;
        TS_te = [TS_te,tets_edge];
        TS_time = mdsvalue(['dim_of(' node ':te_rz)']);
        TS_z = mdsvalue([node ':z_sorted']);
        zts_edge = mdsvalue('\fiber_z');
        TS_z = [TS_z;zts_edge]; 
        if length(zts_edge)~=size(tets_edge,2)
            Te_PF = NaN(size(timebase));
            return
        end
        mdsclose;
        
        Te_PF=NaN(length(TS_time),1);
        itimes = find(TS_time > 0 & TS_time<timebase(end));
         bminor=interp1(efittime,bminor,TS_time, 'linear');
         z0=interp1(efittime,z0,TS_time, 'linear'); 
        for i = 1:length(itimes)
            Te_array=TS_te(itimes(i),:);
            indx=find(Te_array>0);
            if length(indx)<10
                continue
            end
            Te_array=Te_array(indx);
            TS_z_array=TS_z(indx);
            [TS_z_array,index]=sort(TS_z_array);
            Te_array=Te_array(index);
            z_array=linspace(z0(itimes(i)),TS_z_array(end),length(TS_z_array));
            Te_array=interp1(TS_z_array,Te_array,z_array,'linear');
            coreindex = find(z_array<(z0(itimes(i))+0.2*bminor(itimes(i)))&z_array>(z0(itimes(i))-0.2*bminor(itimes(i))));
            if length(coreindex)<2
                continue
            end
            Te_PF(itimes(i))=mean(Te_array(coreindex))/mean(Te_array);
        end
        Te_PF = interp1(TS_time,Te_PF,timebase_column, 'linear');
    else
        Te_PF = NaN(size(timebase));
    end
    
    return
    end
    
    calib=nan;
    mdsopen('electrons',shot);
    TS_ne = mdsvalue([node ':ne_rz']);
    TS_time = mdsvalue(['dim_of(' node ':ne_rz)']);
    TS_z = mdsvalue([node ':z_sorted']);
    nets_edge = mdsvalue('\ts_ne');
    zts_edge = mdsvalue('\fiber_z');
    TS_ne = [TS_ne,nets_edge];
    TS_z = [TS_z;zts_edge];
    if mod(statD,2) == 1&& mod(statG,2) == 1
        TS_te = mdsvalue([node ':te_rz'])*1000*11600;
        tets_edge = mdsvalue('\ts_te')*11600;
        TS_te = [TS_te,tets_edge];
        TS_pressure = TS_te.*TS_ne*1.38e-23;
    end
    mdsclose;
    
    if (mean(nlts1)~=1e32) || (mean(nlts2)~=1e32)
        
        if mean(nlts1)~=1e32 && mean(nlts2)~=1e32
            calib=mean([nltci1;nltci2])/mean([nlts1;nlts2]);
        elseif mean(nlts1)~=1e32
            calib=mean(nltci1)/mean(nlts1);
        else 
            calib=mean(nltci2)/mean(nlts2);
        end
    end
   
    if (calib>0.5) && (calib<1.5)
        TS_ne = TS_ne*calib;
    else
        Te_PF = NaN(size(timebase));
        ne_PF = NaN(size(timebase));
        pressure_PF = NaN(size(timebase));
        return
    end
    
    ne_PF = NaN(length(TS_time),1);
    pressure_PF = NaN(length(TS_time),1);
    Te_PF=NaN(length(TS_time),1);
    itimes = find(TS_time > 0 & TS_time<timebase(end));
    
    if  mod(statD,2)&& mod(statG,2) == 1
         bminor=interp1(efittime,bminor,TS_time, 'linear');
         z0=interp1(efittime,z0,TS_time, 'linear'); 
        for i = 1:length(itimes)
            ne_array=TS_ne(itimes(i),:);
            Te_array=TS_te(itimes(i),:);
            pressure_array=TS_pressure(itimes(i),:);
            indx=find(ne_array>0);
            if length(indx)<10
                continue
            end
            ne_array=ne_array(indx);
            Te_array=Te_array(indx);
            pressure_array=pressure_array(indx);
            TS_z_array=TS_z(indx);
            [TS_z_array,index]=sort(TS_z_array);
            Te_array=Te_array(index);
            ne_array=ne_array(index);
            pressure_array=pressure_array(index);
            z_array=linspace(z0(itimes(i)),TS_z_array(end),length(TS_z_array));
            Te_array=interp1(TS_z_array,Te_array,z_array,'linear');
            ne_array=interp1(TS_z_array,ne_array,z_array,'linear');
            pressure_array=interp1(TS_z_array,pressure_array,z_array,'linear');
            coreindex = find(z_array<(z0(itimes(i))+0.2*bminor(itimes(i)))&z_array>(z0(itimes(i))-0.2*bminor(itimes(i))));
            if length(coreindex)<2
                continue
            end
            ne_PF(itimes(i)) = mean(ne_array(coreindex))/mean(ne_array);
            pressure_PF(itimes(i)) = mean(pressure_array(coreindex))/mean(pressure_array);
            Te_PF(itimes(i))=mean(Te_array(coreindex))/mean(Te_array);
        end
    else
         bminor=interp1(efittime,bminor,TS_time, 'linear');
         z0=interp1(efittime,z0,TS_time, 'linear'); 
        for i = 1:length(itimes)
            ne_array=TS_ne(itimes(i),:);
            indx=find(ne_array>0);
            if length(indx)<10
                continue
            end
            ne_array=ne_array(indx);
            TS_z_array=TS_z(indx);
            [TS_z_array,index]=sort(TS_z_array);
            ne_array=ne_array(index);
            z_array=linspace(z0(itimes(i)),TS_z_array(end),length(TS_z_array));
            ne_array=interp1(TS_z_array,ne_array,z_array,'linear');
            coreindex = find(z_array<(z0(itimes(i))+0.2*bminor(itimes(i)))&z_array>(z0(itimes(i))-0.2*bminor(itimes(i))));
            
            if length(coreindex)<2
                continue
            end
            ne_PF(itimes(i)) = mean(ne_array(coreindex))/mean(ne_array);
        end
    end
    
    ne_PF = interp1(TS_time,ne_PF,timebase_column, 'linear'); 
    pressure_PF = interp1(TS_time, pressure_PF, timebase_column, 'linear');
    Te_PF = interp1(TS_time,Te_PF,timebase_column, 'linear'); 
else
  Te_PF = NaN(size(timebase));
  ne_PF = NaN(size(timebase));
  pressure_PF = NaN(size(timebase));
  return;
end


end