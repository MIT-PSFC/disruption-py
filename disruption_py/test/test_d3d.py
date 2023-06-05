import pytest

from disruption_py.shots import D3DShot, D3D_DISRUPTED_SHOT
from disruption_py.database import create_d3d_handler 

@pytest.fixture 
def d3d():
    return create_d3d_handler()

@pytest.fixture
def d3d_shot():
    shot = D3DShot(D3D_DISRUPTED_SHOT, 'EFIT05',
                disruption_time=4.369214483261109)

def test_init_shot(d3d_shot):
    assert isinstance(d3d_shot, D3DShot)

def test_timebase(d3d_shot):
    assert len(d3d_shot._times) == 2368