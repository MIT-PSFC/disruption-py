"""
Resources
---------
https://docs.python.org/3/library/profile.html

"""

import profile
import pstats
import argparse
from pstats import SortKey

from disruption_py.shots import D3DShot, CModShot, get_shot_id_type
from disruption_py.database import create_d3d_handler, create_cmod_handler


def calibrate_profiler(pr):
    return pr.calibrate(1000000)


def main(args):
    shot_type = get_shot_id_type(args.shot)
    if shot_type == CModShot:
        cmod = create_cmod_handler()
        loop = lambda : CModShot(args.shot)
    elif shot_type == D3DShot:
        d3d = create_d3d_handler()
        loop = lambda : D3DShot(args.shot, args.efit)
    else:
        raise ValueError(f"Shot {args.shot} not found in CMOD or D3D")
    for i in range(args.num_calls):
        loop()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Profile a Shot')
    parser.add_argument('--shot', type=int, default=175552)
    parser.add_argument('--efit', type=str, default='EFIT05')
    parser.add_argument('--num_calls', type=int, default=1)
    args = parser.parse_args()
    pr = profile.Profile()
    pr.bias = calibrate_profiler(pr)
    pr.run('main(args)')
    p = pstats.Stats(pr)
    p.sort_stats(SortKey.CUMULATIVE).print_stats(20)
    p.sort_stats(SortKey.TIME).print_stats(20)
    pr.dump_stats(f"shot_{args.shot}_{args.efit}_profileX{args.num_calls}")
