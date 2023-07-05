"""
Resources
---------
https://docs.python.org/3/library/profile.html

"""

import profile
import pstats
import argparse
import time
from pstats import SortKey

from disruption_py.shots import D3DShot, CModShot, get_shot_id_type
from disruption_py.database import create_d3d_handler, create_cmod_handler
from disruption_py.shots import D3DShot, CModShot, get_shot_id_type
from disruption_py.database import create_d3d_handler, create_cmod_handler

TEST_SHOTS = [1150805012,   # Flattop Disruption
            1150805013,     # No Disruption
            1150805014,     # No Disruption
            1150805015,     # Rampdown Disruption
            1150805016,     # Rampdown Disruption
            1150805017,     # Rampdown Disruption
            1150805019,     # Rampdown Disruption
            1150805020,     # Rampdown Disruption
            1150805021,     # Rampdown Disruption
            1150805022]     # Flattop Disruption

def calibrate_profiler(pr):
    return pr.calibrate(1000000)

def run_cmod_perf():
    cmod = create_cmod_handler()
    times = []
    start_time = time.time()
    for shot_id in TEST_SHOTS:
        shot = CModShot(shot_id, disruption_time=cmod.get_disruption_time(shot_id))
        times.append(time.time() - start_time)
        start_time = time.time()
    return times

def main(args):
    if args.cmod_perf:
        times = run_cmod_perf()
        print(times)
    else:
        shot_type = get_shot_id_type(args.shot)
        pr = profile.Profile()
        pr.bias = calibrate_profiler(pr)
        if shot_type == CModShot:
            cmod = create_cmod_handler()
            pr.run('CModShot(args.shot)')
        elif shot_type == D3DShot:
            d3d = create_d3d_handler()
            pr.run('D3DShot(args.shot, args.efit)')
        else:
            raise ValueError(f"Shot {args.shot} not found in CMOD or D3D")
        p = pstats.Stats(pr)
        p.sort_stats(SortKey.CUMULATIVE).print_stats(50)
        p.sort_stats(SortKey.TIME).print_stats(50)
        pr.dump_stats(f"shot_{args.shot}_profileX{args.num_calls}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Profile a Shot')
    parser = argparse.ArgumentParser(description='Profile a Shot')
    parser.add_argument('--shot', type=int, default=175552)
    parser.add_argument('--efit', type=str, default='EFIT05')
    parser.add_argument('--num_calls', type=int, default=1)
    parser.add_argument('--cmod_perf', action='store_true')
    args = parser.parse_args()
    # ch = logging.StreamHandler(sys.stdout)
    # ch.setLevel(args.log_level*10)
    # # log_format = '%(asctime)s | %(levelname)s: %(message)s'
    # # ch.setFormatter(logging.Formatter(log_format))
    # LOGGER.addHandler(ch)
    # LOGGER.setLevel(args.log_level*10)
    main(args)
