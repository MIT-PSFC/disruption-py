import time

import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt

from disruption_py.inout.mds import ProcessMDSConnection
from disruption_py.machine.tokamak import resolve_tokamak_from_environment

tokamak = resolve_tokamak_from_environment("cmod")

NUM_ACTIVE_WIRE_SEGMENTS = 2
Z_WIRE_INDEX = 2


def z_factor(mds, num_mds_requests):
    for i in range(num_mds_requests):
        mds.get_data(
            rf"\dpcs::top.seg_{i % NUM_ACTIVE_WIRE_SEGMENTS + 1:02d}:p_{Z_WIRE_INDEX:02d}:predictor:factor",
            tree_name="hybrid",
        )


def z_factor_MANY(mds, num_mds_requests):
    names = []
    expressions = []
    for i in range(num_mds_requests):
        names.append(str(i))
        expressions.append(
            rf"\dpcs::top.seg_{i % NUM_ACTIVE_WIRE_SEGMENTS + 1:02d}:p_{Z_WIRE_INDEX:02d}:predictor:factor"
        )
    mds.get_many(names=names, expressions=expressions, tree_name="hybrid")


def get_run_times(method, num_trials, num_mds_requests):
    mds = ProcessMDSConnection.from_config(tokamak).get_shot_connection(1150805012)

    times = []
    for _ in range(num_trials):
        start = time.time()
        method(mds, num_mds_requests)
        end = time.time()
        times.append(end - start)
    return times


NUM_TRIALS = 10

many_avgs = []
orig_avgs = []
num_mds_requests = np.arange(10, 300, 10)
for num_req in tqdm(num_mds_requests):
    orig = get_run_times(z_factor, NUM_TRIALS, num_req)
    many = get_run_times(z_factor_MANY, NUM_TRIALS, num_req)

    many = np.array(many) * 1000
    orig = np.array(orig) * 1000

    many_avgs.append(np.mean(many))
    orig_avgs.append(np.mean(orig))

plt.title(f"MDS getMany vs get for different number of z_factors, {NUM_TRIALS} trials")
plt.scatter(num_mds_requests, many_avgs, label="getMany")
plt.scatter(num_mds_requests, orig_avgs, label="get")
plt.xlabel("Number of MDSplus requests")
plt.ylabel("Time (ms)")
plt.legend()
plt.savefig("z_factor_many.png")
plt.show()
