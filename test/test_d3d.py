import pickle

import pytest
import numpy as np

from disruption_py.shots import D3DShot, D3D_DISRUPTED_SHOT

@pytest.fixture(scope="module")
def d3d_shot():
    shot = D3DShot(D3D_DISRUPTED_SHOT, 'EFIT05',
                disruption_time=4.369214483261109)
    return shot 

@pytest.fixture(scope="module")
def test_shots():
    # Open pickle file test_shot_data.pkl
    with open('data/test_shot_data.pkl', 'rb') as f:
        shots = pickle.load(f)
    return shots

def test_init_shot(d3d_shot):
    assert isinstance(d3d_shot, D3DShot)

#TODO: Compare arrays instead of length 
def test_default_timebase(d3d_shot):
    assert len(d3d_shot._times) == 2368

def test_disruption_timebase():
    shot = D3DShot(D3D_DISRUPTED_SHOT, 'EFIT05',
                disruption_time=4.369214483261109, timebase_signal='disruption_timebase')
    disruption_time = np.load('data/test_shot_disruption_times.npy')
    np.testing.assert_allclose(shot._times, disruption_time)

def test_random_timebase():
    shot = D3DShot(D3D_DISRUPTED_SHOT, 'EFIT05',
                disruption_time=4.369214483261109, timebase_signal='ip')
    ip_time = np.load('data/test_shot_times.npy')
    np.testing.assert_allclose(shot._times, ip_time)

# TODO: Add tests for all parameter methods. Can use the the data in test_shot_data.pkl(contains data from SQL for 10 shots) to test. 
def test_ip_parameters(test_shots):
    pass 

def test_h_parameters(test_shots):
    pass 
