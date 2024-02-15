'''
Methods to load and plot C-Mod TS data

maris, May 2023 (based on sciortino, August 2020)
'''
import matplotlib.pyplot as plt
plt.ion()
import numpy as np
import sys
import scipy as sp
import traceback
import logging 
try:
    sys.path.append('/home/sciortino/usr/python3modules/eqtools3')
    sys.path.append('/home/sciortino/usr/python3modules/profiletools3')
    sys.path.append('/home/sciortino/usr/python3modules/gptools3')
    import eqtools
    import profiletools
except Exception as e:
    logging.warning('Could not import profiletools or eqtools')
    logging.debug(traceback.format_exc())
    pass

shot = 1140515017#1090826016#-->LDL shot #114013018 #1140515017
timebase=np.arange(0.5,1.06,0.01)
rhobase=np.arange(0,1,0.001)
edge_rho_min=0.85
edge_rho_max=0.95
# require edge Thomson to be available
p_Te= profiletools.Te(int(shot), include=['CTS','ETS'], abscissa='sqrtpsinorm',
			t_min=timebase.min(),t_max=timebase.max(), remove_zeros=True)
p_ne= profiletools.ne(int(shot), include=['CTS','ETS'], abscissa='sqrtpsinorm',
			t_min=timebase.min(),t_max=timebase.max(), remove_zeros=True)
#p_Te= profiletools.Te(int(shot), include=['ETS'], abscissa='r/a',t_min=tmin,t_max=tmax, remove_zeros=True)
#p_ne= profiletools.ne(int(shot), include=['ETS'], abscissa='r/a',t_min=tmin,t_max=tmax, remove_zeros=True)


#try:
#	equal_R = p_ne.X[:,1] == p_Te.X[:,1]
#	print(p_ne.X[:,1])
#	print(p_Te.X[:,1])
#	print(equal_R)
#	print(len(p_ne.X[:,1]))
#	print(len(p_Te.X[:,1]))
#	assert np.sum(equal_R) == len(p_ne.X[:,1])
#except:
#	raise ValueError('Edge Thomson rhobase differs between ne and Te')


# consider only flux surface on which points were measured, regardless of LFS or HFS
p_Te.X=np.abs(p_Te.X)
p_ne.X=np.abs(p_ne.X)

# set some minimum uncertainties. Recall that units in objects are 1e20m^{-3} and keV
p_ne.y[p_ne.y<=0.] = 0.01  # 10^18 m^-3
p_Te.y[p_Te.y<=0.01] = 0.01 # 10 eV
p_ne.err_y[p_ne.err_y<=0.01] = 0.01 # 10^18 m^-3
p_Te.err_y[p_Te.err_y<=0.02] = 0.02 # 20 eV

# points in the pedestal that have x uncertainties larger than 0.1 don't help at all
# do this filtering here because filtering of err_X only works before time-averaging
p_ne.remove_points(np.logical_and(p_ne.X[:,1]>0.9, p_ne.err_X[:,1]>0.1))
p_Te.remove_points(np.logical_and(p_Te.X[:,1]>0.9, p_Te.err_X[:,1]>0.1))

# cleanup of low Te values
p_Te.remove_points(np.logical_and(p_Te.X[:,0]<1.03, p_Te.y<0.015))  # TS Te should be >15 eV inside near SOL

#Linear interpolate on time and rho
Te_interpolator = sp.interpolate.LinearNDInterpolator((p_Te.X[:,0],p_Te.X[:,1]), p_Te.y) 
ne_interpolator = sp.interpolate.LinearNDInterpolator((p_ne.X[:,0],p_ne.X[:,1]), p_ne.y) 
#Create mesh to compute interpolation over
timebase_mesh, rhobase_mesh = np.meshgrid(timebase,rhobase)
#Compute interpolated values
Te_interp = Te_interpolator(timebase_mesh,rhobase_mesh)
ne_interp = ne_interpolator(timebase_mesh,rhobase_mesh)

plotting=True
if plotting:
	plt.pcolormesh(timebase_mesh, rhobase_mesh, Te_interp, shading='auto')
	plt.plot(p_Te.X[:,0],p_Te.X[:,1], "ok", label="input point,Te")
	plt.plot(p_ne.X[:,0],p_ne.X[:,1], "ok", c='blue', label="input point,Ne")
	plt.legend()
	plt.colorbar()
	plt.show(block=True)

if plotting:
	plt.pcolormesh(timebase_mesh, rhobase_mesh, ne_interp, shading='auto')
	plt.plot(p_Te.X[:,0],p_Te.X[:,1], "ok", label="input point,Te")
	plt.plot(p_ne.X[:,0],p_ne.X[:,1], "ok", c='blue', label="input point,Ne")
	plt.legend()
	plt.colorbar()
	plt.show(block=True)

if plotting:
	plt.pcolormesh(timebase_mesh, rhobase_mesh, ne_interp, shading='auto')
	#plt.plot(p_Te.X[:,0],p_Te.X[:,1], "ok", label="input point,Te")
	#plt.plot(p_ne.X[:,0],p_ne.X[:,1], "ok", c='blue', label="input point,Ne")
	plt.legend()
	plt.colorbar()
	plt.show(block=True)

###########Compute Te_edge
#Make mask for rho in edge region
rhobase_mesh_mask = (rhobase_mesh >= edge_rho_min) & (rhobase_mesh <= edge_rho_max)#((rhobase_mesh >= edge_rho_min) & (rhobase_mesh <= edge_rho_max))
#Assert that rho values are indeed in desired range
assert np.all(rhobase_mesh[rhobase_mesh_mask]>= edge_rho_min)
assert np.all(rhobase_mesh[rhobase_mesh_mask]<= edge_rho_max)
#Use mask to get only edge values

#Replace all non-edge values with nans
Te_interp_edge=np.where(rhobase_mesh_mask,Te_interp,np.nan) 
ne_interp_edge=np.where(rhobase_mesh_mask,ne_interp,np.nan)
#ne_interp_edge=ne_interp[rhobase_mesh_mask].reshape((num_edge_points,len(timebase)))
#Compute edge quantities
#print(f'rhobase_mesh.shape {rhobase_mesh.shape}')
#print(f'rhobase_mesh_mask {rhobase_mesh_mask}')
#print(f'rhobase_mesh_mask {rhobase_mesh_mask}')
#print(f'Te_interp.shape, {Te_interp.shape}')
#print(f'Te_interp_edge.shape. {Te_interp_edge.shape}')

if plotting:
	plt.pcolormesh(timebase_mesh, rhobase_mesh, Te_interp_edge, shading='auto')
	plt.plot(p_Te.X[:,0],p_Te.X[:,1], "ok", label="input point,Te")
	plt.legend()
	plt.colorbar()
	plt.show(block=True)



Te_edge = np.nanmean(Te_interp_edge,axis=0) 
ne_edge = np.nanmean(ne_interp_edge,axis=0) 
assert len(Te_edge) == len(timebase)
assert len(ne_edge) == len(timebase)

if plotting:
	plt.plot(timebase,Te_edge,label='Te_edge')
	plt.plot(timebase,ne_edge,label='ne_edge')
	plt.legend()
	plt.show(block=True)




#lam_T_nl = 1 # placeholder
## use two point model to get T_sep
##Te_sep_eV = fit_2D.Teu_2pt_model(shot,tmin,tmax,np.mean(in_lambdaT), geqdsk, pressure_opt = 3, lambdaq_opt=1)
#Te_sep_eV = fit_2D.Teu_2pt_model(shot, tmin, tmax, lam_T_nl, geqdsk, pressure_opt = 3, lambdaq_opt=1)
#print('Brunner Te LCFS eV', Te_sep_eV)


#[...]


#nu_star = calc_nu_star(ne, Te)










