import math 

import pandas as pd
import numpy as np 
from scipy.interpolate import interp1d

from MDSplus import *

DEFAULT_SHOT_COLUMNS = ['time','shot','time_until_disrupt','ip']

class Shot(dict):
    """
    Base shot class to represent a single shot of a fusion experiment.

    Attributes
    ---
    metadata->dict: 
        labels: {} /* Dictionary where key is the label name(e.g. "disruption period")
                    and value is a list of starting and ending time tuples(e.g.[(10,100)]*/
        commit_hash: // Commit hash of code version that was used to generate shot data 
        timestep: // Time delta. Every array must conform to this time step
        duration: // Length in milliseconds of shot
        description: // Optional text summary to add to shot
    data->pandas.DataFrame:
    """

    # TODO: Populate metadata dict
    def __init__(self, mdsplus_name,shot_id,data_columns = DEFAULT_SHOT_COLUMNS,data = None):
        self._tree = Tree(mdsplus_name,shot_id)
        self._shot_id = shot_id
        self._metadata = {
            'labels': {},
            'commit_hash': {},
            'timestep': {},
            'duration': {},
            'description': "",
            'disrupted': 100  #TODO: Fix
        }
        try:
            self._analysis_tree = Tree('efit18', shot_id, mode="readonly")
            self._times = self._efit18_tree.getNode(r"\efit18::efit_aeqdsk:time").getData().data()
        except mdsExceptions.TreeFOPENR as e:
            try:
                self._analysis_tree = Tree('analysis', shot_id, mode="readonly")
                self._times = self._analysis_tree.getNode(r"\analysis::efit_aeqdsk:time").getData().data()
            except mdsExceptions.TreeFOPENR as f: 
                print("WARNING: No EFIT18 data found")
        self.data = pd.DataFrame()
        self._populate_shot_data()
    
    def _populate_shot_data(self):
        self.data['times'] = self._times
        self.data['time_until_disrupt'] = self._calc_time_until_disrupt() 
        self.data['ip'],self.data['dip'] = self._calc_Ip_parameters()
        self.data['kappa_area'] = self._calc_kappa_area()
    
    def _calc_time_until_disrupt(self):
        if self._metadata['disrupted']:
            if self._metadata['disrupted'] == -1:
                raise ValueError("Shot disrupted but no t_disrupt found.")   
            time_until_disrupt = self._metadata['disrupted'] - self.data['times'] 
        return time_until_disrupt

    # TODO: Fix plasma current calculation. Add programmed current
    def _calc_Ip_parameters(self):
        magnetics_tree = Tree('magnetics',self._shot_id)
        ip = magnetics_tree.getNode(r"\ip").getData().data()
        magtime = magnetics_tree.getNode(r"\ip").getData().getDims()[0]
        print(magtime.__dict__)
        dip = np.gradient(ip, magtime)
        set_ip_interp = interp1d(magtime, ip, kind='linear')
        set_dip_interp = interp1d(magtime, dip,kind='linear')
        ip = set_ip_interp(self._times)
        dip = set_dip_interp(self._times)
        return ip, dip

    def _calc_Z_parameters(self):
        pass

    def _calc_power(self):
        pass

    def _calc_EFIT_parameters(self):
        pass
    
    def _calc_kappa_area(self):
        aminor = self._analysis_tree.getNode(r'\efit_a_eqdsk:aminor').getData().data()
        area = self._analysis_tree.getNode(r'\efit_a_eqdsk:area').getData().data()
        kappa_area = area/(np.pi * aminor**2)
        set_kappa_interp = interp1d(self._times, kappa_area, kind='linear')
        kappa_area = set_kappa_interp(self._times) #Redundant for now I know, this is a placeholder
        return kappa_area

    def _calc_rotation_velocity(self):
        pass
    
    def _calc_n_equal_1_amplitude(self):
        pass
         
    def _calc_TS_data_cmod(self):
        pass

    def _calc_densities(self):
        pass

    def _calc_efc_current(self):
        pass

    def _calc_sxr_data(self):
        pass        

    def __getitem__(self, key):
        return self._metadata if key == 'metadata' else self.data[key]

class D3DShot(Shot):
    pass

class EASTShot(Shot):
    pass

class CMODShot(Shot):
    pass 
    