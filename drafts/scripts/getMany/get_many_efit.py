import time

from matplotlib import pyplot as plt
import numpy as np
from tqdm import tqdm
from disruption_py.settings.retrieval_settings import RetrievalSettings
from disruption_py.workflow import get_shots_data
from tests.test_cmod_writeback import assert_frame_equal_unordered


shots = [1150805012, 1150805013, 1150805014]


def get_efit():
    return get_shots_data(
        tokamak="cmod",
        shotlist_setting=shots,
        retrieval_settings=RetrievalSettings(
            run_methods=["_get_EFIT_parameters"],
            run_tags=[],
        ),
        output_setting="dataframe",
        num_processes=1,
    )


def get_many_efit():
    return get_shots_data(
        tokamak="cmod",
        shotlist_setting=shots,
        retrieval_settings=RetrievalSettings(
            run_methods=["_getMany_EFIT_parameters"],
            run_tags=[],
        ),
        output_setting="dataframe",
        num_processes=1,
    )


def get_run_times(method, num_trials):
    times = []
    for _ in tqdm(range(num_trials)):
        start = time.time()
        method()
        end = time.time()
        times.append(end - start)
    return times


result = get_many_efit()
assert_frame_equal_unordered(get_many_efit(), get_efit())


NUM_TRIALS = 100

many_avgs = []
orig_avgs = []
orig = get_run_times(get_efit, NUM_TRIALS)
many = get_run_times(get_many_efit, NUM_TRIALS)

orig = np.array(orig) * 1000
many = np.array(many) * 1000

plt.title(f"MDS getMany vs get for efit requests, {NUM_TRIALS} trials, 22 quantities")
plt.hist(orig, label="get", alpha=0.7)
plt.hist(many, label="getMany", alpha=0.7)

orig_mean, many_mean = np.mean(orig), np.mean(many)
plt.axvline(orig_mean, linestyle="--", color="C0", label=f"$\mu={orig_mean: .2f}$ ms")
plt.axvline(many_mean, linestyle="--", color="C1", label=f"$\mu={many_mean: .2f}$ ms")
plt.xlabel("Time to get efit data (ms)")
plt.ylabel("Frequency")
plt.legend()
plt.tight_layout()
plt.savefig("efit_many.png")
