from unittest.mock import patch
import pytest

from disruption_py.utils.mappings.tokamak_helpers import get_tokamak_from_environment, get_tokamak_test_expected_failure_columns, get_tokamak_handler, get_tokamak_test_shot_ids, get_tokamak_test_columns
from disruption_py.utils.math_utils import matlab_gradient_1d_vectorized

def pytest_addoption(parser):
    parser.addoption("--verbose_output", action="store_true", help="More testing information.")
    parser.addoption("--fail_quick", action="store_true", help="Finish test and report statistics instead of failing fast.")

@pytest.fixture(scope="session")
def verbose_output(pytestconfig):
    return pytestconfig.getoption("verbose_output")

@pytest.fixture(scope="session")
def fail_quick(pytestconfig):
    return pytestconfig.getoption("fail_quick")

def pytest_generate_tests(metafunc):
    tokamak = get_tokamak_from_environment()
    
    # parameterized across tests
    if "data_column" in metafunc.fixturenames:
        test_columns = get_tokamak_test_columns(tokamak)
        metafunc.parametrize("data_column", test_columns)
    
@pytest.fixture(scope="session")
def tokamak():
    return get_tokamak_from_environment()

@pytest.fixture(scope="module")
def handler(tokamak):
    return get_tokamak_handler(tokamak)

@pytest.fixture(scope="module")
def shotlist(tokamak):
    return get_tokamak_test_shot_ids(tokamak)

@pytest.fixture(scope="module")
def data_columns(tokamak):
    return get_tokamak_test_columns(tokamak)

@pytest.fixture(scope="module")
def expected_failure_columns(tokamak):
    return get_tokamak_test_expected_failure_columns(tokamak)

# for testing against sql, values generated with matlab use a different gradient method that must be patched for testing
@pytest.fixture(scope='session', autouse=True)
def mock_numpy_gradient():
    with patch('numpy.gradient', new=matlab_gradient_1d_vectorized):
        # The patch will be in place for the duration of the test session
        yield