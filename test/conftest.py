from unittest.mock import patch
import pytest

from disruption_py.utils.math_utils import matlab_gradient_1d_vectorized

@pytest.fixture(scope='session', autouse=True)
def mock_numpy_gradient():
    with patch('numpy.gradient', new=matlab_gradient_1d_vectorized):
        # The patch will be in place for the duration of the test session
        yield


def pytest_addoption(parser):
    parser.addoption("--verbose_output", action="store_true", help="More testing information.")
    parser.addoption("--fail_slow", action="store_true", help="Finish test and report statistics instead of failing fast.")

@pytest.fixture(scope="session")
def verbose_output(pytestconfig):
    return pytestconfig.getoption("verbose_output")

@pytest.fixture(scope="session")
def fail_slow(pytestconfig):
    return pytestconfig.getoption("fail_slow")