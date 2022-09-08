function [ne_PF, Te_PF, pressure_PF] = get_peaking_factor_cmod(shot, timebase)


%path('/home/jxzhu/matlab/yags', path);
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
% [~, status] = mdsopen('cmod',shot);  % Already opened in calling program
status = 1;
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
%      mdsclose;
       return
   end
% mdsclose;
end

%[~, status] = mdsopen('electrons',shot);
[nnodes, status] = mdsvalue('getnci("\\top.electrons","number_of_children")');
if (mod(status,2)==0 || nnodes == 0);
    ne_PF = NaN(size(timebase));
    pressure_PF = NaN(size(timebase));
    Te_PF = NaN(size(timebase));
%   mdsclose;
    return
end
node = '\top.electrons.yag_new.results.profiles';

[~,statA] = mdsvalue([node ':ne_rz']);
[~,statB] = mdsvalue(['dim_of(' node ':ne_rz)']);
[~,statC] = mdsvalue([node ':z_sorted']);
[~,statD] = mdsvalue([node ':te_rz']);
nlnum=4;

if (mod(statA,2) == 1 && mod(statB,2) == 1 && mod(statC,2) == 1)
    try
    [nlts1, nlts2, nltci1, nltci2,~,~]=compare_ts_tci(shot, nlnum);
    catch
        ne_PF = NaN(size(timebase));
        pressure_PF = NaN(size(timebase));
%       mdsopen('electrons',shot);
        if mod(statD,2) == 1
        TS_te = mdsvalue([node ':te_rz'])*1000*1.6022e-19;
        TS_time = mdsvalue(['dim_of(' node ':te_rz)']);
        TS_z = mdsvalue([node ':z_sorted']);
        end
%       mdsclose;
        
    
    if a && mod(statD,2) == 1
        Te_PF=NaN(length(TS_time),1);
        itimes = find(TS_time > 0);
         bminor=interp1(efittime,bminor,TS_time, 'linear');
         z0=interp1(efittime,z0,TS_time, 'linear'); 
        for i = 1:length(itimes)
            coreindex = TS_z<(z0(itimes(i))+0.15*bminor(itimes(i)))&TS_z>(z0(itimes(i))-0.15*bminor(itimes(i)));
            Te_PF(itimes(i))=mean(TS_te(itimes(i),coreindex))/mean(TS_te(itimes(i),:));
        end
        Te_PF = interp1(TS_time,Te_PF,timebase_column, 'linear');
    else
        Te_PF = NaN(size(timebase));
    end
    
    return
    end
    calib=nan;
%   mdsopen('electrons',shot);
    TS_ne = mdsvalue([node ':ne_rz']);
    TS_time = mdsvalue(['dim_of(' node ':ne_rz)']);
    TS_z = mdsvalue([node ':z_sorted']);
    if mod(statD,2) == 1
        TS_te = mdsvalue([node ':te_rz'])*1000*1.6022e-19;
        TS_pressure = TS_te.*TS_ne;
    end
%   mdsclose;
    
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
    itimes = find(TS_time > 0);
    
    if a && mod(statD,2) == 1
         bminor=interp1(efittime,bminor,TS_time, 'linear');
         z0=interp1(efittime,z0,TS_time, 'linear'); 
        for i = 1:length(itimes)
            coreindex = TS_z<(z0(itimes(i))+0.15*bminor(itimes(i)))&TS_z>(z0(itimes(i))-0.15*bminor(itimes(i)));
            ne_PF(itimes(i)) = mean(TS_ne(itimes(i),coreindex))/mean(TS_ne(itimes(i),:));
            pressure_PF(itimes(i)) = mean(TS_pressure(itimes(i),coreindex))/mean(TS_pressure(itimes(i),:));
            Te_PF(itimes(i))=mean(TS_te(itimes(i),coreindex))/mean(TS_te(itimes(i),:));
        end
    elseif a == 1
         bminor=interp1(efittime,bminor,TS_time, 'linear');
         z0=interp1(efittime,z0,TS_time, 'linear'); 
        for i = 1:length(itimes)
            coreindex = TS_z<(z0(itimes(i))+0.15*bminor(itimes(i)))&TS_z>(z0(itimes(i))-0.15*bminor(itimes(i)));
            ne_PF(itimes(i)) = mean(TS_ne(itimes(i),coreindex))/mean(TS_ne(itimes(i),:));
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