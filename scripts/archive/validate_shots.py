import argparse
import random

from disruption_py.database import *
from disruption_py.shots import *

TEST_D3D_SHOT_LIST = [175552, 175553, 175554, 191914, "191786"]

if __name__ == "__main__":
    logger = logging.getLogger("disruption_py")
    # Output to terminal
    # ch = logging.StreamHandler()

    # Output to file:
    ch = logging.FileHandler("./output/validation.log")

    ch.setLevel(logging.DEBUG)
    logger.addHandler(ch)

    parser = argparse.ArgumentParser(
        description="Validate shots against their information in the disruption databases."
    )
    parser.add_argument(
        "shot_list",
        nargs="*",
        help="List of shots to validate. If no shots are provided, all shots in the database will be validated.",
    )
    parser.add_argument(
        "--num_shots",
        type=int,
        help="Specify a number of random shots to validate on.",
        default=None,
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        help="Output directory for plots.",
        default="./output/",
    )
    parser.add_argument(
        "--visualize", type=bool, help="Show comparison plots", default=False
    )
    args = parser.parse_args()

    d3d = create_d3d_handler()
    shot_list = TEST_D3D_SHOT_LIST
    if args.shot_list:
        shot_list = args.shot_list
    elif args.num_shots is not None:
        shot_list = random.choices(
            d3d.get_disruption_shotlist()["shot"].values.tolist(), k=args.num_shots
        )
        print(shot_list)
    validation_results = []
    for shot in shot_list:
        validation_results.append(
            d3d.validate_shot(
                shot, visualize_differences=args.visualize, output_dir=args.output_dir
            )
        )
    for result, shot in zip(validation_results, shot_list):
        print(f"Shot {shot} validation result: {result[0]}")
        if isinstance(result[1], str):
            print(result[1])
