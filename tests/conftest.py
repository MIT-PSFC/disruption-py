from unittest.mock import patch
import pytest

from disruption_py.utils.eval.environment_constants import get_test_expected_failure_columns, get_test_handler, get_test_shot_ids, get_test_columns
from disruption_py.utils.mappings.tokamak_helpers import get_tokamak_from_environment
from disruption_py.utils.math_utils import matlab_gradient_1d_vectorized


def pytest_addoption(parser):
    parser.addoption("--verbose_output", action="store_true", help="More testing information.")
    parser.addoption("--fail_quick", action="store_true", help="Finish test and report statistics instead of failing fast.")

def pytest_generate_tests(metafunc):
    tokamak = get_tokamak_from_environment()
    
    # constants
    if "handler" in metafunc.fixturenames:
        test_handler = get_test_handler(tokamak)    
        metafunc.parametrize("handler", [test_handler], scope="module")
    
    if "shotlist" in metafunc.fixturenames:
        test_shotlist = get_test_shot_ids(tokamak)
        metafunc.parametrize("shotlist", [test_shotlist], scope="module")
    
    # parameterized across tests
    if "data_column" in metafunc.fixturenames:
        test_columns = get_test_columns(tokamak)
        metafunc.parametrize("data_column", test_columns)
        
    if "data_columns" in metafunc.fixturenames:
        test_columns = get_test_columns(tokamak)
        metafunc.parametrize("data_columns", [test_columns])
        
    if "expected_failure_columns" in metafunc.fixturenames:
        test_expected_failure_columns = get_test_expected_failure_columns(tokamak)
        metafunc.parametrize("expected_failure_columns", [test_expected_failure_columns], scope="module")


@pytest.fixture(scope='session', autouse=True)
def mock_numpy_gradient():
    with patch('numpy.gradient', new=matlab_gradient_1d_vectorized):
        # The patch will be in place for the duration of the test session
        yield

@pytest.fixture(scope="session")
def verbose_output(pytestconfig):
    return pytestconfig.getoption("verbose_output")

@pytest.fixture(scope="session")
def fail_quick(pytestconfig):
    return pytestconfig.getoption("fail_quick")