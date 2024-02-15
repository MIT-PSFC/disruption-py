from enum import Enum

'''
For documentation of supported tokamaks:
# --8<-- [start:allowed_tokamak_types_snippet]
Currently supported tokamak type strings are: `"cmod"`
# --8<-- [end:allowed_tokamak_types_snippet]
'''
    
class Tokamak(Enum):
    D3D = 'd3d'
    CMOD = 'cmod'
    EAST = 'east'    
