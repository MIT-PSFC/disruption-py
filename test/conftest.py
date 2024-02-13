import numpy as np
from unittest.mock import patch
import pytest

def matlab_gradient_1d_vectorized(f, h, **kwargs):
    """
    Compute the gradient for a 1D array using vectorized operations.

    :param f: Input 1D array
    :param h: Spacing array with the same length as f
    :return: Gradient of f
    """
    f = np.array(f)
    h = np.array(h)
    
    h_diff = np.diff(h)
    f_diff = np.diff(f)
    # Combine into a single gradient array
    g = np.empty_like(f)
    g[0] = f_diff[0] / h_diff[0]
    g[-1] = f_diff[-1] / h_diff[-1]
    g[1:-1] = (f[2:] - f[0:-2]) / (h[2:] - h[0:-2])

    return g

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