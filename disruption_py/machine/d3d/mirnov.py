#!/usr/bin/env python3

"""
Module for retrieving and processing EFIT parameters for CMOD.
"""

import numpy as np
from MDSplus import mdsExceptions
from scipy.signal import ShortTimeFFT
from scipy.signal.windows import kaiser
import xarray as xr

from disruption_py.core.physics_method.decorator import physics_method, cache_method
from disruption_py.core.physics_method.params import PhysicsMethodParams
from disruption_py.core.utils.math import interp1
from disruption_py.machine.tokamak import Tokamak
from disruption_py.machine.cmod.mirnov import setup_fft

# DIII-D Magnetic Sensors Dictionary
# Data extracted from "Magnetic Diagnostics – Coordinates Of The Sensors" document
# Ted Strait and S. Munaretto - Last Update 2020-05-05
D3D_PROBES_BP = {
    # Bp - 322 degree poloidal array
    'MPI11M322': {'R': 0.973, 'Z': -0.002, 'phi': 322.5, 'gamma': 89.9, 'L': 0.115, 'W': 0.0556, 'NA': 0.04977, 'type': 'Bp'},
    'MPI1A322': {'R': 0.974, 'Z': 0.182, 'phi': 322.5, 'gamma': 90.0, 'L': 0.140, 'W': 0.0556, 'NA': 0.05798, 'type': 'Bp'},
    'MPI2A322': {'R': 0.974, 'Z': 0.512, 'phi': 322.5, 'gamma': 89.7, 'L': 0.140, 'W': 0.0556, 'NA': 0.05929, 'type': 'Bp'},
    'MPI3A322': {'R': 0.975, 'Z': 0.850, 'phi': 322.5, 'gamma': 90.2, 'L': 0.140, 'W': 0.0556, 'NA': 0.05843, 'type': 'Bp'},
    'MPI4A322': {'R': 0.972, 'Z': 1.161, 'phi': 322.5, 'gamma': 90.5, 'L': 0.140, 'W': 0.0556, 'NA': 0.05851, 'type': 'Bp'},
    'MPI5A322': {'R': 1.051, 'Z': 1.330, 'phi': 322.5, 'gamma': 44.6, 'L': 0.141, 'W': 0.0556, 'NA': 0.05891, 'type': 'Bp'},
    'MPI8A322': {'R': 1.219, 'Z': 1.406, 'phi': 322.5, 'gamma': 0.3, 'L': 0.139, 'W': 0.0556, 'NA': 0.06072, 'type': 'Bp'},
    'MPI89A322': {'R': 1.402, 'Z': 1.407, 'phi': 322.5, 'gamma': 0.7, 'L': 0.115, 'W': 0.0556, 'NA': 0.04833, 'type': 'Bp'},
    'MPI9A322': {'R': 1.584, 'Z': 1.408, 'phi': 322.5, 'gamma': 0.6, 'L': 0.140, 'W': 0.0556, 'NA': 0.05917, 'type': 'Bp'},
    'MPI79FA322': {'R': 1.783, 'Z': 1.323, 'phi': 322.5, 'gamma': -39.3, 'L': 0.141, 'W': 0.0556, 'NA': 0.05963, 'type': 'Bp'},
    'MPI79NA322': {'R': 1.924, 'Z': 1.206, 'phi': 322.5, 'gamma': -39.2, 'L': 0.114, 'W': 0.0556, 'NA': 0.04792, 'type': 'Bp'},
    'MPI7FA322': {'R': 2.067, 'Z': 1.090, 'phi': 322.5, 'gamma': -39.4, 'L': 0.140, 'W': 0.0556, 'NA': 0.05700, 'type': 'Bp'},
    'MPI7NA322': {'R': 2.219, 'Z': 0.870, 'phi': 323.3, 'gamma': -68.1, 'L': 0.140, 'W': 0.0556, 'NA': 0.05863, 'type': 'Bp'},
    'MPI67A322': {'R': 2.270, 'Z': 0.746, 'phi': 321.7, 'gamma': -67.9, 'L': 0.141, 'W': 0.0556, 'NA': 0.05998, 'type': 'Bp'},
    'MPI6FA322': {'R': 2.319, 'Z': 0.623, 'phi': 323.2, 'gamma': -68.0, 'L': 0.140, 'W': 0.0556, 'NA': 0.05768, 'type': 'Bp'},
    'MPI6NA322': {'R': 2.416, 'Z': 0.249, 'phi': 317.4, 'gamma': -89.3, 'L': 0.141, 'W': 0.0556, 'NA': 0.05958, 'type': 'Bp'},
    'MPI66M322': {'R': 2.418, 'Z': -0.001, 'phi': 317.4, 'gamma': -89.9, 'L': 0.140, 'W': 0.0556, 'NA': 0.06040, 'type': 'Bp'},
    'MPI1B322': {'R': 0.974, 'Z': -0.187, 'phi': 322.5, 'gamma': 90.1, 'L': 0.140, 'W': 0.0556, 'NA': 0.06017, 'type': 'Bp'},
    'MPI2B322': {'R': 0.975, 'Z': -0.512, 'phi': 322.5, 'gamma': 90.3, 'L': 0.140, 'W': 0.0556, 'NA': 0.05962, 'type': 'Bp'},
    'MPI3B322': {'R': 0.974, 'Z': -0.854, 'phi': 322.5, 'gamma': 89.8, 'L': 0.140, 'W': 0.0556, 'NA': 0.05960, 'type': 'Bp'},
    'MPI4B322': {'R': 0.972, 'Z': -1.159, 'phi': 322.5, 'gamma': 89.6, 'L': 0.141, 'W': 0.0556, 'NA': 0.05856, 'type': 'Bp'},
    'MPI5B322': {'R': 1.048, 'Z': -1.330, 'phi': 322.5, 'gamma': 136.4, 'L': 0.140, 'W': 0.0556, 'NA': 0.06048, 'type': 'Bp'},
    'MPI8B322': {'R': 1.254, 'Z': -1.405, 'phi': 322.5, 'gamma': -180.0, 'L': 0.140, 'W': 0.0556, 'NA': 0.05967, 'type': 'Bp'},
    'MPI89B322': {'R': 1.477, 'Z': -1.406, 'phi': 322.5, 'gamma': -179.9, 'L': 0.140, 'W': 0.0556, 'NA': 0.05817, 'type': 'Bp'},
    'MPI9B322': {'R': 1.699, 'Z': -1.406, 'phi': 322.5, 'gamma': 179.9, 'L': 0.142, 'W': 0.0556, 'NA': 0.05849, 'type': 'Bp'},
    'MPI79B322': {'R': 1.894, 'Z': -1.333, 'phi': 322.5, 'gamma': -129.5, 'L': 0.141, 'W': 0.0556, 'NA': 0.05751, 'type': 'Bp'},
    'MPI7FB322': {'R': 2.085, 'Z': -1.102, 'phi': 322.5, 'gamma': -129.2, 'L': 0.140, 'W': 0.0556, 'NA': 0.05773, 'type': 'Bp'},
    'MPI7NB322': {'R': 2.212, 'Z': -0.873, 'phi': 323.4, 'gamma': -113.5, 'L': 0.141, 'W': 0.0556, 'NA': 0.05741, 'type': 'Bp'},
    'MPI67B322': {'R': 2.263, 'Z': -0.749, 'phi': 321.7, 'gamma': -112.8, 'L': 0.140, 'W': 0.0556, 'NA': 0.06145, 'type': 'Bp'},
    'MPI6FB322': {'R': 2.315, 'Z': -0.624, 'phi': 323.3, 'gamma': -113.4, 'L': 0.140, 'W': 0.0556, 'NA': 0.05904, 'type': 'Bp'},
    'MPI6NB322': {'R': 2.416, 'Z': -0.244, 'phi': 317.4, 'gamma': -90.8, 'L': 0.140, 'W': 0.0556, 'NA': 0.06041, 'type': 'Bp'},

    # Bp - 67 degree backup probes
    'MPI2A067': {'R': 0.973, 'Z': 0.518, 'phi': 75.4, 'gamma': 89.8, 'L': 0.140, 'W': 0.0556, 'NA': 0.06015, 'type': 'Bp'},
    'MPI11M067': {'R': 0.973, 'Z': -0.004, 'phi': 75.4, 'gamma': 90.0, 'L': 0.141, 'W': 0.0556, 'NA': 0.06095, 'type': 'Bp'},
    'MPI2B067': {'R': 0.972, 'Z': -0.517, 'phi': 75.4, 'gamma': 90.3, 'L': 0.139, 'W': 0.0556, 'NA': 0.05997, 'type': 'Bp'},
    'MPI67A097': {'R': 2.266, 'Z': 0.758, 'phi': 97.4, 'gamma': -67.6, 'L': 0.155, 'W': 0.1000, 'NA': 0.13076, 'type': 'Bp'},
    'MPI66M067': {'R': 2.413, 'Z': 0.003, 'phi': 67.5, 'gamma': -89.9, 'L': 0.141, 'W': 0.0556, 'NA': 0.06218, 'type': 'Bp'},
    'MPI67B097': {'R': 2.262, 'Z': -0.751, 'phi': 97.4, 'gamma': -112.6, 'L': 0.156, 'W': 0.1000, 'NA': 0.13311, 'type': 'Bp'},

    # Bp - 142 degree poloidal array
    'MPI1A139': {'R': 0.979, 'Z': 0.070, 'phi': 142.5, 'gamma': 90.3, 'L': 0.142, 'W': 0.0556, 'NA': 0.09776, 'type': 'Bp'},
    'MPI2A139': {'R': 0.977, 'Z': 0.209, 'phi': 135.0, 'gamma': 89.9, 'L': 0.141, 'W': 0.0556, 'NA': 0.09808, 'type': 'Bp'},
    'MPI3A139': {'R': 0.979, 'Z': 0.347, 'phi': 142.5, 'gamma': 90.1, 'L': 0.141, 'W': 0.0556, 'NA': 0.09725, 'type': 'Bp'},
    'MPI4A139': {'R': 0.976, 'Z': 0.485, 'phi': 135.0, 'gamma': 90.0, 'L': 0.141, 'W': 0.0556, 'NA': 0.09767, 'type': 'Bp'},
    'MPI5A139': {'R': 0.979, 'Z': 0.759, 'phi': 142.5, 'gamma': 89.7, 'L': 0.141, 'W': 0.0556, 'NA': 0.09756, 'type': 'Bp'},
    'MPI79A147': {'R': 1.762, 'Z': 1.312, 'phi': 147.2, 'gamma': -39.5, 'L': 0.153, 'W': 0.1000, 'NA': 0.21266, 'type': 'Bp'},
    'MPI67A142': {'R': 2.262, 'Z': 0.757, 'phi': 142.7, 'gamma': -67.3, 'L': 0.153, 'W': 0.1000, 'NA': 0.13527, 'type': 'Bp'},
    'MPI67A157': {'R': 2.264, 'Z': 0.760, 'phi': 157.5, 'gamma': -67.8, 'L': 0.156, 'W': 0.1000, 'NA': 0.13329, 'type': 'Bp'},
    'MPI6NA132': {'R': 2.409, 'Z': 0.277, 'phi': 132.7, 'gamma': -89.0, 'L': 0.138, 'W': 0.0556, 'NA': 0.06170, 'type': 'Bp'},
    'MPI6NA157': {'R': 2.407, 'Z': 0.254, 'phi': 157.6, 'gamma': -89.4, 'L': 0.141, 'W': 0.0556, 'NA': 0.05951, 'type': 'Bp'},
    'MPI66M157': {'R': 2.413, 'Z': -0.001, 'phi': 157.6, 'gamma': -89.8, 'L': 0.137, 'W': 0.0556, 'NA': 0.06109, 'type': 'Bp'},
    'MPI6NB157': {'R': 2.414, 'Z': -0.244, 'phi': 157.6, 'gamma': -89.6, 'L': 0.141, 'W': 0.0556, 'NA': 0.06336, 'type': 'Bp'},
    'MPI6FB142': {'R': 2.312, 'Z': -0.620, 'phi': 143.4, 'gamma': -113.0, 'L': 0.140, 'W': 0.0556, 'NA': 0.06147, 'type': 'Bp'},
    'MPI67B157': {'R': 2.260, 'Z': -0.754, 'phi': 157.4, 'gamma': -112.8, 'L': 0.155, 'W': 0.1000, 'NA': 0.13229, 'type': 'Bp'},
    'MPI7NB142': {'R': 2.209, 'Z': -0.867, 'phi': 143.3, 'gamma': -112.6, 'L': 0.142, 'W': 0.0556, 'NA': 0.06023, 'type': 'Bp'},
    'MPI79B142': {'R': 2.044, 'Z': -1.110, 'phi': 142.4, 'gamma': -129.1, 'L': 0.154, 'W': 0.1000, 'NA': 0.20995, 'type': 'Bp'},
    'MPI5B139': {'R': 0.977, 'Z': -0.757, 'phi': 134.9, 'gamma': 89.8, 'L': 0.142, 'W': 0.0556, 'NA': 0.09822, 'type': 'Bp'},
    'MPI4B139': {'R': 0.980, 'Z': -0.482, 'phi': 142.5, 'gamma': 89.9, 'L': 0.141, 'W': 0.0556, 'NA': 0.09748, 'type': 'Bp'},
    'MPI3B139': {'R': 0.977, 'Z': -0.343, 'phi': 135.0, 'gamma': 90.1, 'L': 0.141, 'W': 0.0556, 'NA': 0.09816, 'type': 'Bp'},
    'MPI2B139': {'R': 0.979, 'Z': -0.207, 'phi': 142.5, 'gamma': 89.6, 'L': 0.141, 'W': 0.0556, 'NA': 0.09812, 'type': 'Bp'},
    'MPI1B139': {'R': 0.977, 'Z': -0.069, 'phi': 135.0, 'gamma': 90.0, 'L': 0.142, 'W': 0.0556, 'NA': 0.09748, 'type': 'Bp'},
    'MPI1B157': {'R': 0.972, 'Z': -0.183, 'phi': 157.5, 'gamma': 89.6, 'L': 0.141, 'W': 0.0556, 'NA': 0.06068, 'type': 'Bp'},

    # Bp - upper divertor baffles
    'MPI1U157': {'R': 1.445, 'Z': 1.277, 'phi': 155.0, 'gamma': -40.6, 'L': 0.027, 'W': 0.0556, 'NA': 0.01198, 'type': 'Bp'},
    'MPI2U157': {'R': 1.560, 'Z': 1.187, 'phi': 155.0, 'gamma': -40.8, 'L': 0.027, 'W': 0.0556, 'NA': 0.01227, 'type': 'Bp'},
    'MPI3U157': {'R': 1.723, 'Z': 1.112, 'phi': 155.0, 'gamma': -3.3, 'L': 0.027, 'W': 0.0556, 'NA': 0.01229, 'type': 'Bp'},
    'MPI4U157': {'R': 1.873, 'Z': 1.116, 'phi': 155.4, 'gamma': -2.8, 'L': 0.054, 'W': 0.0556, 'NA': 0.02398, 'type': 'Bp'},
    'MPI5U157': {'R': 1.218, 'Z': 1.298, 'phi': 152.2, 'gamma': 66.3, 'L': 0.027, 'W': 0.0556, 'NA': 0.01139, 'type': 'Bp'},
    'MPI6U157': {'R': 1.075, 'Z': 1.202, 'phi': 152.2, 'gamma': 0.7, 'L': 0.027, 'W': 0.0556, 'NA': 0.01154, 'type': 'Bp'},
    'MPI7U157': {'R': 0.976, 'Z': 0.979, 'phi': 157.5, 'gamma': 90.3, 'L': 0.027, 'W': 0.0556, 'NA': 0.01144, 'type': 'Bp'},

    # Bp - R0 ELM array (127-137 deg.)
    'MPI66M127': {'R': 2.413, 'Z': 0.002, 'phi': 127.9, 'gamma': -90.0, 'L': 0.140, 'W': 0.0556, 'NA': 0.05984, 'type': 'Bp'},
    'MPI66M132': {'R': 2.409, 'Z': 0.007, 'phi': 132.5, 'gamma': -90.3, 'L': 0.055, 'W': 0.0556, 'NA': 0.02387, 'type': 'Bp'},
    'MPI66M137': {'R': 2.416, 'Z': 0.006, 'phi': 137.4, 'gamma': -90.4, 'L': 0.054, 'W': 0.0556, 'NA': 0.02310, 'type': 'Bp'},
    'MPI66B137': {'R': 2.410, 'Z': -0.124, 'phi': 137.5, 'gamma': -89.8, 'L': 0.053, 'W': 0.0556, 'NA': 0.02448, 'type': 'Bp'},
    'MPI6NB137': {'R': 2.409, 'Z': -0.255, 'phi': 137.5, 'gamma': -90.5, 'L': 0.053, 'W': 0.0556, 'NA': 0.02378, 'type': 'Bp'},

    # Bp - R0 ELM array (307-317 deg)
    'MPI66M307': {'R': 2.413, 'Z': 0.001, 'phi': 307.0, 'gamma': -90.2, 'L': 0.140, 'W': 0.0556, 'NA': 0.05984, 'type': 'Bp'},
    'MPI66M312': {'R': 2.411, 'Z': -0.001, 'phi': 312.4, 'gamma': -90.0, 'L': 0.054, 'W': 0.0556, 'NA': 0.02297, 'type': 'Bp'},
    'MPI6NA312': {'R': 2.411, 'Z': 0.264, 'phi': 313.1, 'gamma': -89.2, 'L': 0.055, 'W': 0.0556, 'NA': 0.02419, 'type': 'Bp'},
    'MPI66B312': {'R': 2.411, 'Z': -0.131, 'phi': 313.2, 'gamma': -90.3, 'L': 0.055, 'W': 0.0556, 'NA': 0.02398, 'type': 'Bp'},
    'MPI6NB312': {'R': 2.411, 'Z': -0.261, 'phi': 313.2, 'gamma': -90.3, 'L': 0.055, 'W': 0.0556, 'NA': 0.02394, 'type': 'Bp'},

    # Bp - lower divertor toroidal array
    'MPI1L020': {'R': 1.5070, 'Z': -1.2885, 'phi': 20.0, 'gamma': -178.17, 'L': 0.0271, 'W': 0.0556, 'NA': 0.01207, 'type': 'Bp'},
    'MPI2L020': {'R': 1.6064, 'Z': -1.2892, 'phi': 20.0, 'gamma': -179.69, 'L': 0.0274, 'W': 0.0556, 'NA': 0.01212, 'type': 'Bp'},
    'MPI1L050': {'R': 1.5055, 'Z': -1.2883, 'phi': 50.0, 'gamma': -179.73, 'L': 0.0269, 'W': 0.0556, 'NA': 0.01198, 'type': 'Bp'},
    'MPI1L110': {'R': 1.5127, 'Z': -1.2900, 'phi': 110.0, 'gamma': -178.96, 'L': 0.0279, 'W': 0.0556, 'NA': 0.01226, 'type': 'Bp'},
    'MPI1L180': {'R': 1.5104, 'Z': -1.2915, 'phi': 180.0, 'gamma': -179.99, 'L': 0.0271, 'W': 0.0556, 'NA': 0.01130, 'type': 'Bp'},
    'MPI2L180': {'R': 1.6103, 'Z': -1.2915, 'phi': 180.0, 'gamma': -179.75, 'L': 0.0277, 'W': 0.0556, 'NA': 0.01117, 'type': 'Bp'},
    'MPI3L180': {'R': 1.7103, 'Z': -1.2918, 'phi': 180.0, 'gamma': -178.98, 'L': 0.0273, 'W': 0.0556, 'NA': 0.01125, 'type': 'Bp'},
    'MPI1L230': {'R': 1.5070, 'Z': -1.2920, 'phi': 230.0, 'gamma': 179.61, 'L': 0.0271, 'W': 0.0556, 'NA': 0.01103, 'type': 'Bp'},
    'MPI1L320': {'R': 1.5055, 'Z': -1.2899, 'phi': 320.0, 'gamma': 179.74, 'L': 0.0273, 'W': 0.0556, 'NA': 0.01194, 'type': 'Bp'},

    # Bp - R0 toroidal array
    'MPI66M020': {'R': 2.410, 'Z': 0.002, 'phi': 19.5, 'gamma': -89.9, 'L': 0.139, 'W': 0.0556, 'NA': 0.06202, 'type': 'Bp'},
    'MPI66M097': {'R': 2.413, 'Z': -0.005, 'phi': 97.4, 'gamma': -89.8, 'L': 0.137, 'W': 0.0556, 'NA': 0.06087, 'type': 'Bp'},
    'MPI66M200': {'R': 2.412, 'Z': 0.003, 'phi': 199.7, 'gamma': -89.9, 'L': 0.139, 'W': 0.0556, 'NA': 0.06048, 'type': 'Bp'},
    'MPI66M247': {'R': 2.413, 'Z': -0.003, 'phi': 246.4, 'gamma': -90.5, 'L': 0.140, 'W': 0.0556, 'NA': 0.06150, 'type': 'Bp'},
    'MPI66M277': {'R': 2.413, 'Z': -0.009, 'phi': 277.5, 'gamma': -89.7, 'L': 0.137, 'W': 0.0556, 'NA': 0.06194, 'type': 'Bp'},
    'MPI66M340': {'R': 2.413, 'Z': -0.002, 'phi': 339.7, 'gamma': -89.8, 'L': 0.140, 'W': 0.0556, 'NA': 0.06072, 'type': 'Bp'},

    # Bp - R+1 toroidal array
    'MPI67A022': {'R': 2.273, 'Z': 0.767, 'phi': 22.3, 'gamma': -67.5, 'L': 0.153, 'W': 0.1000, 'NA': 0.14410, 'type': 'Bp'},
    'MPI67A037': {'R': 2.268, 'Z': 0.755, 'phi': 37.5, 'gamma': -67.9, 'L': 0.155, 'W': 0.1000, 'NA': 0.13018, 'type': 'Bp'},
    'MPI67A052': {'R': 2.270, 'Z': 0.764, 'phi': 52.4, 'gamma': -65.9, 'L': 0.154, 'W': 0.1000, 'NA': 0.14439, 'type': 'Bp'},
    'MPI67A067': {'R': 2.261, 'Z': 0.758, 'phi': 67.5, 'gamma': -68.0, 'L': 0.153, 'W': 0.1000, 'NA': 0.13683, 'type': 'Bp'},
    'MPI67A082': {'R': 2.260, 'Z': 0.757, 'phi': 82.5, 'gamma': -67.6, 'L': 0.154, 'W': 0.1000, 'NA': 0.13749, 'type': 'Bp'},
    'MPI67A217': {'R': 2.268, 'Z': 0.753, 'phi': 217.5, 'gamma': -67.6, 'L': 0.155, 'W': 0.1000, 'NA': 0.13148, 'type': 'Bp'},
    'MPI67A262': {'R': 2.264, 'Z': 0.751, 'phi': 262.6, 'gamma': -67.5, 'L': 0.154, 'W': 0.1000, 'NA': 0.13937, 'type': 'Bp'},
    'MPI67A277': {'R': 2.267, 'Z': 0.756, 'phi': 277.3, 'gamma': -68.0, 'L': 0.155, 'W': 0.1000, 'NA': 0.13266, 'type': 'Bp'},
    'MPI67A307': {'R': 2.266, 'Z': 0.760, 'phi': 307.5, 'gamma': -67.5, 'L': 0.152, 'W': 0.1000, 'NA': 0.13519, 'type': 'Bp'},
    'MPI67A337': {'R': 2.268, 'Z': 0.756, 'phi': 337.4, 'gamma': -67.6, 'L': 0.155, 'W': 0.1000, 'NA': 0.13068, 'type': 'Bp'},

    # Bp - R-1 toroidal array
    'MPI67B022': {'R': 2.253, 'Z': -0.751, 'phi': 22.4, 'gamma': -112.3, 'L': 0.154, 'W': 0.1000, 'NA': 0.14287, 'type': 'Bp'},
    'MPI67B037': {'R': 2.264, 'Z': -0.752, 'phi': 37.4, 'gamma': -112.4, 'L': 0.156, 'W': 0.1000, 'NA': 0.13078, 'type': 'Bp'},
    'MPI67B052': {'R': 2.257, 'Z': -0.751, 'phi': 52.4, 'gamma': -112.7, 'L': 0.154, 'W': 0.1000, 'NA': 0.14125, 'type': 'Bp'},
    'MPI67B217': {'R': 2.266, 'Z': -0.753, 'phi': 217.6, 'gamma': -112.8, 'L': 0.155, 'W': 0.1000, 'NA': 0.13554, 'type': 'Bp'},
    'MPI67B277': {'R': 2.260, 'Z': -0.756, 'phi': 277.5, 'gamma': -112.5, 'L': 0.155, 'W': 0.1000, 'NA': 0.13036, 'type': 'Bp'},
    'MPI67B337': {'R': 2.260, 'Z': -0.754, 'phi': 337.4, 'gamma': -112.7, 'L': 0.155, 'W': 0.1000, 'NA': 0.14188, 'type': 'Bp'},

    # Bp - R+2 toroidal array
    'MPI79A072': {'R': 1.762, 'Z': 1.312, 'phi': 72.2, 'gamma': -39.0, 'L': 0.153, 'W': 0.1000, 'NA': 0.21191, 'type': 'Bp'},
    'MPI79A222': {'R': 1.762, 'Z': 1.312, 'phi': 222.1, 'gamma': -39.5, 'L': 0.153, 'W': 0.1000, 'NA': 0.21037, 'type': 'Bp'},
    'MPI79A272': {'R': 1.763, 'Z': 1.311, 'phi': 272.9, 'gamma': -39.1, 'L': 0.153, 'W': 0.1000, 'NA': 0.21067, 'type': 'Bp'},

    # Bp - R-2 toroidal array
    'MPI79B067': {'R': 2.039, 'Z': -1.110, 'phi': 67.5, 'gamma': -129.4, 'L': 0.153, 'W': 0.1000, 'NA': 0.21217, 'type': 'Bp'},
    'MPI79B217': {'R': 2.046, 'Z': -1.109, 'phi': 217.5, 'gamma': -129.6, 'L': 0.153, 'W': 0.1000, 'NA': 0.21082, 'type': 'Bp'},
    'MPI79B277': {'R': 2.043, 'Z': -1.110, 'phi': 277.6, 'gamma': -129.2, 'L': 0.153, 'W': 0.1000, 'NA': 0.21190, 'type': 'Bp'},

    # Bp - HFS vertical array (139 deg)
    'MPI5A139': {'R': 0.979, 'Z': 0.759, 'phi': 142.5, 'gamma': 89.7, 'L': 0.141, 'W': 0.0556, 'NA': 0.09756, 'type': 'Bp'},

    # Bp - HFS vertical array (199 deg)
    'MPI5A199': {'R': 0.978, 'Z': 0.759, 'phi': 202.6, 'gamma': 90.0, 'L': 0.141, 'W': 0.0556, 'NA': 0.09805, 'type': 'Bp'},
    'MPI4A199': {'R': 0.976, 'Z': 0.485, 'phi': 195.0, 'gamma': 90.2, 'L': 0.141, 'W': 0.0556, 'NA': 0.09741, 'type': 'Bp'},
    'MPI3A199': {'R': 0.977, 'Z': 0.347, 'phi': 202.6, 'gamma': 90.0, 'L': 0.141, 'W': 0.0556, 'NA': 0.09803, 'type': 'Bp'},
    'MPI2A199': {'R': 0.976, 'Z': 0.209, 'phi': 195.0, 'gamma': 89.8, 'L': 0.141, 'W': 0.0556, 'NA': 0.09762, 'type': 'Bp'},
    'MPI1A199': {'R': 0.976, 'Z': 0.070, 'phi': 202.6, 'gamma': 90.6, 'L': 0.141, 'W': 0.0556, 'NA': 0.09796, 'type': 'Bp'},
    'MPI1B199': {'R': 0.976, 'Z': -0.069, 'phi': 195.0, 'gamma': 90.3, 'L': 0.141, 'W': 0.0556, 'NA': 0.09842, 'type': 'Bp'},
    'MPI2B199': {'R': 0.977, 'Z': -0.207, 'phi': 202.6, 'gamma': 89.8, 'L': 0.141, 'W': 0.0556, 'NA': 0.09740, 'type': 'Bp'},
    'MPI3B199': {'R': 0.976, 'Z': -0.344, 'phi': 195.0, 'gamma': 90.0, 'L': 0.142, 'W': 0.0556, 'NA': 0.09782, 'type': 'Bp'},
    'MPI4B199': {'R': 0.978, 'Z': -0.483, 'phi': 202.5, 'gamma': 89.7, 'L': 0.142, 'W': 0.0556, 'NA': 0.09836, 'type': 'Bp'},
    'MPI5B199': {'R': 0.977, 'Z': -0.757, 'phi': 194.9, 'gamma': 90.0, 'L': 0.142, 'W': 0.0556, 'NA': 0.09831, 'type': 'Bp'},

    # Bp - HFS toroidal array (Above midplane)
    'MPI1A011': {'R': 0.977, 'Z': 0.071, 'phi': 15.1, 'gamma': 90.3, 'L': 0.141, 'W': 0.0556, 'NA': 0.09765, 'type': 'Bp'},
    'MPI1A049': {'R': 0.979, 'Z': 0.070, 'phi': 52.6, 'gamma': 90.2, 'L': 0.141, 'W': 0.0556, 'NA': 0.09753, 'type': 'Bp'},
    'MPI1A109': {'R': 0.977, 'Z': 0.070, 'phi': 112.5, 'gamma': 90.4, 'L': 0.141, 'W': 0.0556, 'NA': 0.09834, 'type': 'Bp'},
    'MPI1A244': {'R': 0.977, 'Z': 0.070, 'phi': 247.3, 'gamma': 89.9, 'L': 0.141, 'W': 0.0556, 'NA': 0.09751, 'type': 'Bp'},
    'MPI1A274': {'R': 0.978, 'Z': 0.071, 'phi': 277.5, 'gamma': 90.3, 'L': 0.141, 'W': 0.0556, 'NA': 0.09771, 'type': 'Bp'},
    'MPI1A341': {'R': 0.978, 'Z': 0.071, 'phi': 345.1, 'gamma': 89.9, 'L': 0.141, 'W': 0.0556, 'NA': 0.09756, 'type': 'Bp'},

    # Bp - HFS toroidal array (Below midplane)
    'MPI1B011': {'R': 0.976, 'Z': -0.068, 'phi': 7.6, 'gamma': 90.1, 'L': 0.141, 'W': 0.0556, 'NA': 0.09785, 'type': 'Bp'},
    'MPI1B049': {'R': 0.978, 'Z': -0.068, 'phi': 45.1, 'gamma': 90.0, 'L': 0.141, 'W': 0.0556, 'NA': 0.09860, 'type': 'Bp'},
    'MPI1B109': {'R': 0.976, 'Z': -0.069, 'phi': 105.0, 'gamma': 89.8, 'L': 0.141, 'W': 0.0556, 'NA': 0.09854, 'type': 'Bp'},
    'MPI1B244': {'R': 0.977, 'Z': -0.069, 'phi': 239.9, 'gamma': 90.3, 'L': 0.141, 'W': 0.0556, 'NA': 0.09777, 'type': 'Bp'},
    'MPI1B274': {'R': 0.976, 'Z': -0.070, 'phi': 270.1, 'gamma': 90.2, 'L': 0.141, 'W': 0.0556, 'NA': 0.09762, 'type': 'Bp'},
    'MPI1B341': {'R': 0.977, 'Z': -0.069, 'phi': 337.6, 'gamma': 89.8, 'L': 0.141, 'W': 0.0556, 'NA': 0.09755, 'type': 'Bp'},
}

D3D_PROBES_BR = {
    # Br - upper divertor baffles
    'DSL1U180': {'R': 1.506, 'Z': 1.270, 'phi': 185.0, 'gamma': 49.5, 'L': 0.094, 'W': 60.1, 'NA': 0.13912, 'type': 'Br'},
    'DSL2U180': {'R': 1.590, 'Z': 1.198, 'phi': 185.0, 'gamma': 49.5, 'L': 0.088, 'W': 60.2, 'NA': 0.14843, 'type': 'Br'},
    'DSL3U180': {'R': 1.730, 'Z': 1.142, 'phi': 185.2, 'gamma': 92.3, 'L': 0.106, 'W': 59.2, 'NA': 0.18909, 'type': 'Br'},
    'DSL4U157': {'R': 1.877, 'Z': 1.133, 'phi': 157.6, 'gamma': 90.5, 'L': 0.196, 'W': 6.0, 'NA': 0.07450, 'type': 'Br'},
    'DSL5U157': {'R': 1.185, 'Z': 1.286, 'phi': 168.7, 'gamma': 156.6, 'L': 0.091, 'W': 108.0, 'NA': 0.20466, 'type': 'Br'},
    'DSL6U157': {'R': 1.086, 'Z': 1.228, 'phi': 168.7, 'gamma': 90.4, 'L': 0.107, 'W': 107.6, 'NA': 0.21782, 'type': 'Br'},

    # Br - R0 toroidal array
    'ISL66M017': {'R': 2.431, 'Z': 0.000, 'phi': 17.2, 'gamma': 0.0, 'L': 0.800, 'W': 60.3, 'NA': 2.0466, 'type': 'Br'},
    'ISL66M042': {'R': 2.431, 'Z': 0.000, 'phi': 40.5, 'gamma': 0.0, 'L': 0.800, 'W': 56.0, 'NA': 1.9040, 'type': 'Br'},
    'ISL66M072': {'R': 2.431, 'Z': 0.000, 'phi': 72.4, 'gamma': 0.0, 'L': 0.800, 'W': 50.1, 'NA': 1.7007, 'type': 'Br'},
    'ISL66M102': {'R': 2.431, 'Z': 0.000, 'phi': 100.2, 'gamma': 0.0, 'L': 0.800, 'W': 63.3, 'NA': 2.1522, 'type': 'Br'},
    'ISL66M132': {'R': 2.431, 'Z': 0.000, 'phi': 132.6, 'gamma': 0.0, 'L': 0.800, 'W': 70.3, 'NA': 2.3850, 'type': 'Br'},
    'ISL66M197': {'R': 2.431, 'Z': 0.000, 'phi': 197.6, 'gamma': 0.0, 'L': 0.800, 'W': 59.6, 'NA': 2.0238, 'type': 'Br'},
    'ISL66M252': {'R': 2.431, 'Z': 0.000, 'phi': 252.4, 'gamma': 0.0, 'L': 0.800, 'W': 50.1, 'NA': 1.7007, 'type': 'Br'},
    'ISL66M312': {'R': 2.431, 'Z': 0.000, 'phi': 312.3, 'gamma': 0.0, 'L': 0.800, 'W': 69.6, 'NA': 2.3633, 'type': 'Br'},

    # Br - R+1 toroidal array
    'ISL67A017': {'R': 2.300, 'Z': 0.714, 'phi': 17.8, 'gamma': 22.6, 'L': 0.680, 'W': 61.5, 'NA': 1.6779, 'type': 'Br'},
    'ISL67A052': {'R': 2.300, 'Z': 0.714, 'phi': 50.1, 'gamma': 22.6, 'L': 0.680, 'W': 60.0, 'NA': 1.6371, 'type': 'Br'},
    'ISL67A072': {'R': 2.300, 'Z': 0.714, 'phi': 73.0, 'gamma': 22.6, 'L': 0.680, 'W': 48.9, 'NA': 1.3344, 'type': 'Br'},
    'ISL67A112': {'R': 2.300, 'Z': 0.714, 'phi': 111.3, 'gamma': 22.6, 'L': 0.680, 'W': 62.4, 'NA': 1.7021, 'type': 'Br'},
    'ISL67A132': {'R': 2.300, 'Z': 0.714, 'phi': 133.0, 'gamma': 22.6, 'L': 0.680, 'W': 71.1, 'NA': 1.9396, 'type': 'Br'},
    'ISL67A197': {'R': 2.300, 'Z': 0.714, 'phi': 198.6, 'gamma': 22.6, 'L': 0.680, 'W': 60.0, 'NA': 1.6370, 'type': 'Br'},
    'ISL67A252': {'R': 2.300, 'Z': 0.714, 'phi': 253.0, 'gamma': 22.6, 'L': 0.680, 'W': 48.9, 'NA': 1.3344, 'type': 'Br'},
    'ISL67A312': {'R': 2.300, 'Z': 0.714, 'phi': 312.3, 'gamma': 22.6, 'L': 0.680, 'W': 69.6, 'NA': 1.8997, 'type': 'Br'},

    # Br - R-1 toroidal array
    'ISL67B017': {'R': 2.300, 'Z': -0.714, 'phi': 18.7, 'gamma': -22.6, 'L': 0.680, 'W': 59.7, 'NA': 1.6289, 'type': 'Br'},
    'ISL67B052': {'R': 2.300, 'Z': -0.714, 'phi': 50.1, 'gamma': -22.6, 'L': 0.680, 'W': 60.0, 'NA': 1.6371, 'type': 'Br'},
    'ISL67B072': {'R': 2.300, 'Z': -0.714, 'phi': 73.0, 'gamma': -22.6, 'L': 0.680, 'W': 48.9, 'NA': 1.3344, 'type': 'Br'},
    'ISL67B112': {'R': 2.300, 'Z': -0.714, 'phi': 111.4, 'gamma': -22.6, 'L': 0.680, 'W': 62.5, 'NA': 1.7058, 'type': 'Br'},
    'ISL67B132': {'R': 2.300, 'Z': -0.714, 'phi': 133.0, 'gamma': -22.6, 'L': 0.680, 'W': 71.1, 'NA': 1.9396, 'type': 'Br'},
    'ISL67B197': {'R': 2.300, 'Z': -0.714, 'phi': 198.6, 'gamma': -22.6, 'L': 0.680, 'W': 60.0, 'NA': 1.6370, 'type': 'Br'},
    'ISL67B252': {'R': 2.300, 'Z': -0.714, 'phi': 253.0, 'gamma': -22.6, 'L': 0.680, 'W': 48.9, 'NA': 1.3344, 'type': 'Br'},
    'ISL67B312': {'R': 2.300, 'Z': -0.714, 'phi': 313.2, 'gamma': -22.6, 'L': 0.680, 'W': 71.4, 'NA': 1.9486, 'type': 'Br'},

    # Br - R+2 toroidal array
    'ISL79A072': {'R': 1.773, 'Z': 1.326, 'phi': 72.2, 'gamma': 51.0, 'L': 0.206, 'W': 2.82, 'NA': 0.56651, 'type': 'Br'},
    'ISL79A147': {'R': 1.773, 'Z': 1.327, 'phi': 147.2, 'gamma': 51.0, 'L': 0.206, 'W': 2.82, 'NA': 0.56990, 'type': 'Br'},
    'ISL79A222': {'R': 1.773, 'Z': 1.327, 'phi': 222.1, 'gamma': 51.0, 'L': 0.206, 'W': 2.82, 'NA': 0.55330, 'type': 'Br'},
    'ISL79A272': {'R': 1.774, 'Z': 1.325, 'phi': 272.9, 'gamma': 51.0, 'L': 0.206, 'W': 2.82, 'NA': 0.57139, 'type': 'Br'},

    # Br - R-2 toroidal array
    'ISL79B067': {'R': 2.053, 'Z': -1.122, 'phi': 67.5, 'gamma': -40.1, 'L': 0.206, 'W': 2.43, 'NA': 0.56541, 'type': 'Br'},
    'ISL79B142': {'R': 2.058, 'Z': -1.122, 'phi': 142.4, 'gamma': -40.1, 'L': 0.206, 'W': 2.43, 'NA': 0.56309, 'type': 'Br'},
    'ISL79B217': {'R': 2.060, 'Z': -1.121, 'phi': 217.5, 'gamma': -40.1, 'L': 0.206, 'W': 2.43, 'NA': 0.56656, 'type': 'Br'},
    'ISL79B277': {'R': 2.057, 'Z': -1.122, 'phi': 277.6, 'gamma': -40.1, 'L': 0.206, 'W': 2.43, 'NA': 0.57099, 'type': 'Br'},

    # Br - HFS vertical array (139 deg)
    'ISL5A139': {'R': 0.959, 'Z': 0.759, 'phi': 138.7, 'gamma': 179.9, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL4A139': {'R': 0.959, 'Z': 0.486, 'phi': 138.8, 'gamma': 180.0, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL3A139': {'R': 0.960, 'Z': 0.348, 'phi': 138.8, 'gamma': 180.0, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL2A139': {'R': 0.960, 'Z': 0.209, 'phi': 138.7, 'gamma': 179.7, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL1A139': {'R': 0.961, 'Z': 0.070, 'phi': 138.7, 'gamma': 179.9, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL1B139': {'R': 0.960, 'Z': -0.071, 'phi': 138.7, 'gamma': 180.0, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL2B139': {'R': 0.960, 'Z': -0.208, 'phi': 138.7, 'gamma': 179.9, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL3B139': {'R': 0.959, 'Z': -0.345, 'phi': 138.7, 'gamma': 179.9, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL4B139': {'R': 0.960, 'Z': -0.483, 'phi': 138.7, 'gamma': 179.9, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL5B139': {'R': 0.960, 'Z': -0.762, 'phi': 131.1, 'gamma': 180.0, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},

    # Br - HFS vertical array (199 deg)
    'ISL5A199': {'R': 0.959, 'Z': 0.758, 'phi': 198.8, 'gamma': 179.8, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL4A199': {'R': 0.959, 'Z': 0.486, 'phi': 198.8, 'gamma': 180.2, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL3A199': {'R': 0.959, 'Z': 0.347, 'phi': 198.8, 'gamma': 179.9, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL2A199': {'R': 0.959, 'Z': 0.209, 'phi': 198.7, 'gamma': 179.9, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL1A199': {'R': 0.960, 'Z': 0.070, 'phi': 198.8, 'gamma': 180.0, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL1B199': {'R': 0.959, 'Z': -0.071, 'phi': 198.8, 'gamma': 179.9, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL2B199': {'R': 0.959, 'Z': -0.208, 'phi': 198.8, 'gamma': 180.0, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL3B199': {'R': 0.959, 'Z': -0.345, 'phi': 198.8, 'gamma': 179.9, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL4B199': {'R': 0.959, 'Z': -0.483, 'phi': 198.7, 'gamma': 180.0, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL5B199': {'R': 0.959, 'Z': -0.762, 'phi': 198.7, 'gamma': 180.1, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},

    # Br - HFS toroidal array (Above midplane)
    'ISL1A011': {'R': 0.961, 'Z': 0.070, 'phi': 11.3, 'gamma': 179.9, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL1A049': {'R': 0.960, 'Z': 0.070, 'phi': 48.6, 'gamma': 180.3, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL1A109': {'R': 0.959, 'Z': 0.070, 'phi': 108.7, 'gamma': 180.2, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL1A244': {'R': 0.961, 'Z': 0.070, 'phi': 243.6, 'gamma': 179.4, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL1A274': {'R': 0.960, 'Z': 0.070, 'phi': 273.8, 'gamma': 180.0, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL1A341': {'R': 0.960, 'Z': 0.070, 'phi': 341.3, 'gamma': 179.5, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},

    # Br - HFS toroidal array (Below midplane)
    'ISL1B011': {'R': 0.961, 'Z': -0.069, 'phi': 11.4, 'gamma': 180.5, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL1B049': {'R': 0.960, 'Z': -0.070, 'phi': 48.6, 'gamma': 179.8, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL1B109': {'R': 0.958, 'Z': -0.070, 'phi': 108.7, 'gamma': 179.8, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL1B244': {'R': 0.961, 'Z': -0.069, 'phi': 243.6, 'gamma': 180.6, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL1B274': {'R': 0.960, 'Z': -0.070, 'phi': 273.7, 'gamma': 180.2, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},
    'ISL1B341': {'R': 0.960, 'Z': -0.069, 'phi': 341.3, 'gamma': 180.3, 'L': 0.116, 'W': 20.9, 'NA': 0.24262, 'type': 'Br'},

    # Br - 67 degree poloidal array (outer surface of vacuum vessel)
    'DSL12A067': {'R': 0.925, 'Z': 0.343, 'phi': 67.0, 'gamma': 180.0, 'L': 0.682, 'W': 30.0, 'NA': 0.33011, 'type': 'Br'},
    'DSL34A067': {'R': 0.942, 'Z': 1.003, 'phi': 67.0, 'gamma': 176.9, 'L': 0.631, 'W': 30.0, 'NA': 0.31130, 'type': 'Br'},
    'DSL59A067': {'R': 1.320, 'Z': 1.393, 'phi': 67.0, 'gamma': 101.4, 'L': 0.729, 'W': 29.7, 'NA': 0.49906, 'type': 'Br'},
    'DSL79A067': {'R': 1.946, 'Z': 1.263, 'phi': 67.0, 'gamma': 51.0, 'L': 0.639, 'W': 29.6, 'NA': 0.64262, 'type': 'Br'},
    'DSL67A067': {'R': 2.315, 'Z': 0.780, 'phi': 67.0, 'gamma': 22.4, 'L': 0.557, 'W': 29.7, 'NA': 0.66875, 'type': 'Br'},
    'DSL66M052': {'R': 2.433, 'Z': 0.000, 'phi': 52.0, 'gamma': -0.1, 'L': 1.036, 'W': 24.8, 'NA': 1.09124, 'type': 'Br'},
    'DSL67B067': {'R': 2.313, 'Z': -0.781, 'phi': 67.0, 'gamma': -22.4, 'L': 0.530, 'W': 29.9, 'NA': 0.64055, 'type': 'Br'},
    'DSL79B067': {'R': 2.025, 'Z': -1.262, 'phi': 67.0, 'gamma': -39.9, 'L': 0.519, 'W': 29.9, 'NA': 0.54883, 'type': 'Br'},
    'DSL59B067': {'R': 1.405, 'Z': -1.392, 'phi': 67.0, 'gamma': -99.1, 'L': 0.896, 'W': 29.9, 'NA': 0.65748, 'type': 'Br'},
    'DSL34B067': {'R': 0.942, 'Z': -1.003, 'phi': 67.0, 'gamma': -176.9, 'L': 0.631, 'W': 30.0, 'NA': 0.31130, 'type': 'Br'},
    'DSL12B067': {'R': 0.925, 'Z': -0.343, 'phi': 67.0, 'gamma': 180.0, 'L': 0.682, 'W': 30.0, 'NA': 0.33011, 'type': 'Br'},

    # Br - 157 degree poloidal array (outer surface of vacuum vessel)
    'DSL12A157': {'R': 0.925, 'Z': 0.343, 'phi': 157.0, 'gamma': 180.0, 'L': 0.682, 'W': 30.0, 'NA': 0.33011, 'type': 'Br'},
    'DSL34A157': {'R': 0.942, 'Z': 1.003, 'phi': 157.0, 'gamma': 176.9, 'L': 0.631, 'W': 30.0, 'NA': 0.31130, 'type': 'Br'},
    'DSL59A157': {'R': 1.320, 'Z': 1.393, 'phi': 157.0, 'gamma': 101.4, 'L': 0.729, 'W': 30.0, 'NA': 0.50359, 'type': 'Br'},
    'DSL79A157': {'R': 1.944, 'Z': 1.265, 'phi': 157.0, 'gamma': 51.0, 'L': 0.635, 'W': 29.9, 'NA': 0.64307, 'type': 'Br'},
    'DSL67A157': {'R': 2.314, 'Z': 0.783, 'phi': 157.0, 'gamma': 22.4, 'L': 0.557, 'W': 30.1, 'NA': 0.67718, 'type': 'Br'},
    'DSL66M152': {'R': 2.433, 'Z': 0.001, 'phi': 152.0, 'gamma': 0.0, 'L': 1.037, 'W': 29.1, 'NA': 1.28213, 'type': 'Br'},
    'DSL67B157': {'R': 2.313, 'Z': -0.781, 'phi': 157.0, 'gamma': -22.4, 'L': 0.549, 'W': 30.0, 'NA': 0.66533, 'type': 'Br'},
    'DSL79B157': {'R': 2.025, 'Z': -1.262, 'phi': 157.0, 'gamma': -39.9, 'L': 0.523, 'W': 30.0, 'NA': 0.55500, 'type': 'Br'},
    'DSL59B157': {'R': 1.406, 'Z': -1.392, 'phi': 157.0, 'gamma': -99.1, 'L': 0.897, 'W': 30.0, 'NA': 0.65995, 'type': 'Br'},
    'DSL34B157': {'R': 0.942, 'Z': -1.003, 'phi': 157.0, 'gamma': -176.9, 'L': 0.631, 'W': 30.0, 'NA': 0.31130, 'type': 'Br'},
    'DSL12B157': {'R': 0.925, 'Z': -0.343, 'phi': 157.0, 'gamma': 180.0, 'L': 0.682, 'W': 30.0, 'NA': 0.33011, 'type': 'Br'},

    # Br - R+1 toroidal array (outer surface of vacuum vessel)
    'SL67FA345': {'R': 2.289, 'Z': 0.845, 'phi': 345.0, 'gamma': 22.4, 'L': 0.432, 'W': 21.0, 'NA': 0.36133, 'type': 'Br'},
    'SL67NA345': {'R': 2.395, 'Z': 0.587, 'phi': 345.0, 'gamma': 22.4, 'L': 0.124, 'W': 21.0, 'NA': 0.10858, 'type': 'Br'},

    # Br - R0 toroidal array (outer surface of vacuum vessel)
    'SL66A132': {'R': 2.449, 'Z': 0.261, 'phi': 132.0, 'gamma': 3.5, 'L': 0.519, 'W': 13.1, 'NA': 0.28958, 'type': 'Br'},
    'SL66B132': {'R': 2.447, 'Z': -0.262, 'phi': 132.0, 'gamma': -3.8, 'L': 0.522, 'W': 13.1, 'NA': 0.29140, 'type': 'Br'},
    'SL66A312': {'R': 2.449, 'Z': 0.260, 'phi': 312.0, 'gamma': 3.4, 'L': 0.519, 'W': 12.8, 'NA': 0.28481, 'type': 'Br'},
    'SL66B312': {'R': 2.447, 'Z': -0.261, 'phi': 312.0, 'gamma': -3.7, 'L': 0.522, 'W': 12.8, 'NA': 0.28627, 'type': 'Br'},

    # Br - R-1 toroidal array (outer surface of vacuum vessel)
    'SL67NB015': {'R': 2.393, 'Z': -0.588, 'phi': 15.0, 'gamma': -22.4, 'L': 0.124, 'W': 15.0, 'NA': 0.07765, 'type': 'Br'},
    'SL67FB015': {'R': 2.287, 'Z': -0.845, 'phi': 15.0, 'gamma': -22.4, 'L': 0.430, 'W': 15.0, 'NA': 0.25751, 'type': 'Br'},

    # Br - R0 toroidal array (external)
    'ESL66M019': {'R': 2.477, 'Z': 0.000, 'phi': 19.0, 'gamma': 0.0, 'L': 1.194, 'W': 60.0, 'NA': 3.0960, 'type': 'Br'},
    'ESL66M079': {'R': 2.477, 'Z': 0.000, 'phi': 79.0, 'gamma': 0.0, 'L': 1.194, 'W': 60.0, 'NA': 3.0960, 'type': 'Br'},
    'ESL66M139': {'R': 2.477, 'Z': 0.000, 'phi': 139.0, 'gamma': 0.0, 'L': 1.194, 'W': 60.0, 'NA': 3.0960, 'type': 'Br'},
    'ESL66M199': {'R': 2.477, 'Z': 0.000, 'phi': 199.0, 'gamma': 0.0, 'L': 1.194, 'W': 60.0, 'NA': 3.0960, 'type': 'Br'},
    'ESL66M259': {'R': 2.477, 'Z': 0.000, 'phi': 259.0, 'gamma': 0.0, 'L': 1.194, 'W': 60.0, 'NA': 3.0960, 'type': 'Br'},
    'ESL66M319': {'R': 2.477, 'Z': 0.000, 'phi': 319.0, 'gamma': 0.0, 'L': 1.194, 'W': 60.0, 'NA': 3.0960, 'type': 'Br'},

    # Br - R+1 toroidal array (external)
    'ESL67A004': {'R': 2.389, 'Z': 0.806, 'phi': 4.0, 'gamma': 24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67A034': {'R': 2.389, 'Z': 0.806, 'phi': 34.0, 'gamma': 24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67A064': {'R': 2.389, 'Z': 0.806, 'phi': 64.0, 'gamma': 24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67A094': {'R': 2.389, 'Z': 0.806, 'phi': 94.0, 'gamma': 24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67A124': {'R': 2.389, 'Z': 0.806, 'phi': 124.0, 'gamma': 24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67A154': {'R': 2.389, 'Z': 0.806, 'phi': 154.0, 'gamma': 24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67A184': {'R': 2.389, 'Z': 0.806, 'phi': 184.0, 'gamma': 24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67A214': {'R': 2.389, 'Z': 0.806, 'phi': 214.0, 'gamma': 24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67A244': {'R': 2.389, 'Z': 0.806, 'phi': 244.0, 'gamma': 24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67A274': {'R': 2.389, 'Z': 0.806, 'phi': 274.0, 'gamma': 24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67A304': {'R': 2.389, 'Z': 0.806, 'phi': 304.0, 'gamma': 24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67A334': {'R': 2.389, 'Z': 0.806, 'phi': 334.0, 'gamma': 24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},

    # Br - R-1 toroidal array (external)
    'ESL67B004': {'R': 2.389, 'Z': -0.806, 'phi': 4.0, 'gamma': -24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67B034': {'R': 2.389, 'Z': -0.806, 'phi': 34.0, 'gamma': -24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67B064': {'R': 2.389, 'Z': -0.806, 'phi': 64.0, 'gamma': -24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67B094': {'R': 2.389, 'Z': -0.806, 'phi': 94.0, 'gamma': -24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67B124': {'R': 2.389, 'Z': -0.806, 'phi': 124.0, 'gamma': -24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67B154': {'R': 2.389, 'Z': -0.806, 'phi': 154.0, 'gamma': -24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67B184': {'R': 2.389, 'Z': -0.806, 'phi': 184.0, 'gamma': -24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67B214': {'R': 2.389, 'Z': -0.806, 'phi': 214.0, 'gamma': -24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67B244': {'R': 2.389, 'Z': -0.806, 'phi': 244.0, 'gamma': -24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67B274': {'R': 2.389, 'Z': -0.806, 'phi': 274.0, 'gamma': -24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67B304': {'R': 2.389, 'Z': -0.806, 'phi': 304.0, 'gamma': -24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
    'ESL67B334': {'R': 2.389, 'Z': -0.806, 'phi': 334.0, 'gamma': -24.1, 'L': 0.454, 'W': 30.0, 'NA': 2.2697, 'type': 'Br'},
}



class D3DMirnovMethods:
    """
    Class for retrieving and processing Mirnov-related signals for D3D
    """

    @staticmethod
    def get_mirnov_names_and_locations(params: PhysicsMethodParams, debug=False):
        """Get the names and locations of the Mirnov coils in the ANALYSIS tree.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters for the physics method.

        Returns
        -------
        all_mirnov_names : list[str]
            The names of all the Mirnov coils.
        phi_all : numpy.ndarray
            The toroidal angles of the Mirnov coils
        gamma : np.ndarray
            planar offset of the Mirnov coils
        """
        
    @staticmethod
    def get_mirnov_fft(params: PhysicsMethodParams, mirnov_name: str, freq_resolution: float = 100, max_freq: float = 80e3, suffix="D"):
        """Get and the interpolated fft of a Mirnov coil.
        
        This will work best if the params timebase is uniform.
        """
        # There is a difference between E probes and D probes, D probes tend to be the 500 kHz ones.
        path = f"ptdata('{mirnov_name}{suffix}', {params.shot_id})"
        try:
            mirnov_signal, mirnov_times = params.mds_conn.get_data_with_dims(
                        path=path,
                        tree_name="d3d",
                        astype="float32",
            )
            # DIII-D mirnov times are in ms, convert to seconds
            mirnov_times = mirnov_times / 1000.0
            # Get the sampling frequency of the Mirnov signal
            f_mirnov = 1 / np.mean(np.diff(mirnov_times))
            params.logger.verbose(f"Using Mirnov frequency of {f_mirnov} Hz on {mirnov_name}")
            if not np.isclose(f_mirnov, 500e3, rtol=0.01):
                params.logger.warning(f"[Shot {params.shot_id}] Got Mirnov frequency of {f_mirnov} Hz on {mirnov_name}, expected 500 kHz on DIII-D")

            # Set up the fft taking into account the params timebase
            f_timebase = 1 / np.mean(np.diff(params.times))
            sft = setup_fft(f_mirnov=f_mirnov, f_target=f_timebase, frequency_resolution=freq_resolution)
            f_indices = np.where(sft.f < max_freq)
            freqs = sft.f[f_indices]

            mirnov_fft_full = sft.stft(mirnov_signal)[f_indices]

            # Interpolate the STFFTs to the timebase
            fft_times = (sft.delta_t * np.arange(mirnov_fft_full.shape[1])) + mirnov_times[0]
            mirnov_fft_interp = interp1(fft_times, mirnov_fft_full, params.times)

            # Check if the average difference between frequencies is NOT close to the frequency resolution
            if not np.isclose(np.mean(np.diff(freqs)), freq_resolution, atol=1):
                params.logger.warning(f"[Shot {params.shot_id}] The Mirnov frequency resolution is not {freq_resolution} Hz")
            # Replace the frequencies with nice integers (for the sake of consistency)
            freqs = np.arange(0, max_freq, freq_resolution)

            return mirnov_fft_interp, freqs
        except Exception as e:
            return None, None

    @staticmethod
    @physics_method(
        tokamak=Tokamak.D3D,
    )
    def get_all_mirnov_ffts(params: PhysicsMethodParams):
        """Get all FFTs of the available Mirnov coils for this shot.

        Parameters
        ----------
        params : PhysicsMethodParams
            The parameters for the physics method.

        Returns
        -------
        mirnov_ffts : xarray.DataSet
            The FFT of the Mirnov coils.
            Dimensions are probe, frequency, and time.
            Coordinates are probe, frequency, time, phi, theta, and theta_pol.
        """

        all_mirnov_names, phi_all, theta_all, theta_pol_all = D3DMirnovMethods.get_mirnov_names_and_locations(params, debug=True)

        valid_mirnov_ffts = []
        valid_mirnov_names = []
        valid_mirnov_locations = []
        saved_freqs = None

        for mirnov_name, mirnov_phi, mirnov_theta, mirnov_theta_pol in zip(all_mirnov_names, phi_all, theta_all, theta_pol_all):
            mirnov_fft, freqs = D3DMirnovMethods.get_mirnov_fft(params, mirnov_name)
            if mirnov_fft is not None:
                valid_mirnov_ffts.append(mirnov_fft)
                valid_mirnov_names.append(mirnov_name)
                valid_mirnov_locations.append((mirnov_phi, mirnov_theta, mirnov_theta_pol))

            if saved_freqs is None:
                saved_freqs = freqs

        valid_mirnov_ffts = np.expand_dims(valid_mirnov_ffts, axis=0)  # Add a new axis for the probe dimension

        mirnov_ffts_real = xr.DataArray(
            np.array(valid_mirnov_ffts).real,
            dims=("idx", "probe", "frequency", "time"),
            coords={
                "shot": ("idx", [params.shot_id]),
                "probe": list(range(len(valid_mirnov_locations))),
                "probe_name": ("probe", valid_mirnov_names),
                "frequency": saved_freqs,
                "time": params.times,
                "phi": ("probe", [loc[0] for loc in valid_mirnov_locations]),
                "theta": ("probe", [loc[1] for loc in valid_mirnov_locations]),
                "theta_pol": ("probe", [loc[2] for loc in valid_mirnov_locations]),
            },
        )
        mirnov_ffts_imag = xr.DataArray(
            np.array(valid_mirnov_ffts).imag,
            dims=("idx", "probe", "frequency", "time"),
            coords={
                "shot": ("idx", [params.shot_id]),
                "probe": list(range(len(valid_mirnov_locations))),
                "probe_name": ("probe", valid_mirnov_names),
                "frequency": saved_freqs,
                "time": params.times,
                "phi": ("probe", [loc[0] for loc in valid_mirnov_locations]),
                "theta": ("probe", [loc[1] for loc in valid_mirnov_locations]),
                "theta_pol": ("probe", [loc[2] for loc in valid_mirnov_locations]),
            },
        )
        mirnov_ds_real = mirnov_ffts_real.astype(np.float32).to_dataset(name="mirnov_fft_real")
        mirnov_ds_imag = mirnov_ffts_imag.astype(np.float32).to_dataset(name="mirnov_fft_imag")
        mirnov_ds = xr.merge([mirnov_ds_real, mirnov_ds_imag])

        return mirnov_ds
    
    @staticmethod
    @physics_method(
        tokamak=Tokamak.D3D,
    )
    def get_preferred_mirnov_sxx(params: PhysicsMethodParams):
        """Get the Sxx of a single Mirnov coil in the shot, in order of preference from a few consistently 'good' probes"""

        # From the R0 toroidal array
        preferred_mirnov_names = ['MPI66M020', 'MPI66M097', 'MPI66M200', 'MPI66M247', 'MPI66M277', 'MPI66M340']

        for mirnov_name in preferred_mirnov_names:
            mirnov_fft, freqs = D3DMirnovMethods.get_mirnov_fft(params, mirnov_name)
            if mirnov_fft is not None:
                params.logger.info(f"Using {mirnov_name} for Sxx")
                mirnov_sxx = np.abs(mirnov_fft) ** 2
                return xr.DataArray(
                    mirnov_sxx,
                    dims=("frequency", "time"),
                    coords={
                        "shot": params.shot_id,
                        "frequency": freqs,
                        "time": params.times,
                    },
                    attrs={
                        "probe_name": mirnov_name,
                    },
                ).astype(np.float32).to_dataset(name="mirnov_sxx")