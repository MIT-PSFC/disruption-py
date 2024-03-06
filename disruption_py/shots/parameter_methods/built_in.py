from disruption_py.shots.parameter_methods.cmod.basic_parameter_methods import CModEfitRequests, BasicCmodRequests
from disruption_py.shots.parameter_methods.d3d.basic_parameter_methods import BasicD3DRequests
from disruption_py.shots.parameter_methods.d3d.efit_parameter_methods import D3DEfitRequests
DEFAULT_SHOT_DATA_REQUESTS = [
    CModEfitRequests(), BasicCmodRequests(),
    D3DEfitRequests(), BasicD3DRequests(),
]