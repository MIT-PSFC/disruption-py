import argparse

from disruption_py.shots import *
from disruption_py.database import *


TEST_D3D_SHOT_LIST = [175552,175553,175554, 191914, '191786']

if __name__ == '__main__':
    logging.basicConfig(filename='validation.log', level=logging.DEBUG)
    parser = argparse.ArgumentParser(
        description="Validate shots against their information in the disruption databases.")
    parser.add_argument('shot_list', nargs='*',
                        help='List of shots to validate. If no shots are provided, all shots in the database will be validated.')
    args = parser.parse_args()
    d3d = create_d3d_handler()
    shot_list = TEST_D3D_SHOT_LIST
    if args.shot_list:
        shot_list = args.shot_list
    for shot in shot_list:
        d3d.validate_shot(shot, visualize_differences=False)
