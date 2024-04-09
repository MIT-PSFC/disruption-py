from unittest.mock import patch
import pytest

from disruption_py.handlers.cmod_handler import CModHandler
from disruption_py.handlers.d3d_handler import D3DHandler
from disruption_py.utils.constants import CMOD_TEST_COLUMNS, CMOD_TEST_SHOTS, D3D_TEST_COLUMNS, D3D_TEST_SHOTS
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.mappings.tokamak_helpers import get_tokamak_from_environment
from disruption_py.utils.math_utils import matlab_gradient_1d_vectorized


def pytest_addoption(parser):
    parser.addoption("--verbose_output", action="store_true", help="More testing information.")
    parser.addoption("--fail_slow", action="store_true", help="Finish test and report statistics instead of failing fast.")

def get_test_handler(tokamak : Tokamak):
    if tokamak is Tokamak.CMOD:
        return CModHandler()
    elif tokamak is Tokamak.D3D:
        return D3DHandler()
    else:
        raise ValueError("Tokamak {} not supported for this test".format(tokamak))

def get_test_shotlist(tokamak : Tokamak):
    if tokamak == Tokamak.CMOD:
        return CMOD_TEST_SHOTS
    elif tokamak == Tokamak.D3D:
        return D3D_TEST_SHOTS
    else:
        raise ValueError("Tokamak {} not supported for this test".format(tokamak))
    
def get_test_columns(tokamak : Tokamak):
    if tokamak == Tokamak.CMOD:
        return CMOD_TEST_COLUMNS
    elif tokamak == Tokamak.D3D:
        return D3D_TEST_COLUMNS
    else:
        raise ValueError("Tokamak {} not supported for this test".format(tokamak))

def pytest_generate_tests(metafunc):
    tokamak = get_tokamak_from_environment()
    
    # constants
    test_handler = get_test_handler(tokamak)    
    if "handler" in metafunc.fixturenames:
        metafunc.parametrize("handler", [test_handler], scope="module")
    
    test_shotlist = get_test_shotlist(tokamak)
    if "shotlist" in metafunc.fixturenames:
        metafunc.parametrize("shotlist", [test_shotlist], scope="module")
    
    # parameterized across tests
    test_columns = get_test_columns(tokamak)
    if "data_column" in metafunc.fixturenames:
        metafunc.parametrize("data_column", test_columns)
        
    test_columns = get_test_columns(tokamak)
    if "data_columns" in metafunc.fixturenames:
        metafunc.parametrize("data_columns", [test_columns])


@pytest.fixture(scope='session', autouse=True)
def mock_numpy_gradient():
    with patch('numpy.gradient', new=matlab_gradient_1d_vectorized):
        # The patch will be in place for the duration of the test session
        yield

@pytest.fixture(scope="session")
def verbose_output(pytestconfig):
    return pytestconfig.getoption("verbose_output")

@pytest.fixture(scope="session")
def fail_slow(pytestconfig):
    return pytestconfig.getoption("fail_slow")