from disruption_py.shots import *
from disruption_py.database import *

TEST_D3D_SHOT_LIST = [175552, 191914, '191786']

if __name__ == '__main__':
    d3d = create_d3d_handler()
    for shot in TEST_D3D_SHOT_LIST:
        d3d.validate_shot(shot, visualize_differences=True)
