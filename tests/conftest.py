#!/usr/bin/env python3

"""Pytest configuration module for setting up fixtures."""

import os
import time
from tempfile import mkdtemp
from unittest.mock import patch

import pytest

from disruption_py.core.utils.math import matlab_gradient_1d_vectorized
from disruption_py.machine.tokamak import resolve_tokamak_from_environment
from tests.utils.factory import (
    get_tokamak_test_columns,
    get_tokamak_test_expected_failure_columns,
    get_tokamak_test_shotlist,
)


def pytest_addoption(parser):
    """Add custom command-line options for verbose output to pytest."""
    parser.addoption(
        "--verbose_output", action="store_true", help="More testing information."
    )


@pytest.fixture(scope="session")
def verbose_output(pytestconfig):
    """Fixture to retrieve the verbose output option from pytest configuration."""
    return pytestconfig.getoption("verbose_output")


def pytest_generate_tests(metafunc):
    """
    Parametrize `data_column` and mark expected failure columns. Marked columns
    will xfail on assert False and xpass on assert True.
    """
    tokamak = resolve_tokamak_from_environment()

    # parameterized across tests
    if "data_column" in metafunc.fixturenames:
        test_columns = get_tokamak_test_columns(tokamak)
        xfail_columns = get_tokamak_test_expected_failure_columns(tokamak)
        data_columns = []
        for test_col in test_columns:
            if test_col in xfail_columns:
                data_columns.append(pytest.param(test_col, marks=pytest.mark.xfail))
            else:
                data_columns.append(test_col)
        metafunc.parametrize("data_column", data_columns)


@pytest.fixture(scope="session", name="tokamak")
def tokamak_fixture():
    """Fixture to resolve the tokamak from the environment."""
    return resolve_tokamak_from_environment()


@pytest.fixture(scope="module")
def shotlist(tokamak):
    """Fixture to retrieve the test shotlist for the tokamak."""
    return get_tokamak_test_shotlist(tokamak)


@pytest.fixture(scope="module", name="data_columns")
def data_columns_fixture(tokamak):
    """Fixture to retrieve the test data columns for the tokamak."""
    return get_tokamak_test_columns(tokamak)


@pytest.fixture(scope="module")
def expected_failure_columns(tokamak):
    """Fixture to retrieve the expected failure columns for the tokamak."""
    return get_tokamak_test_expected_failure_columns(tokamak)


@pytest.fixture(scope="session", autouse=True)
def mock_numpy_gradient():
    """
    This fixture patches NumPy's gradient function with a MATLAB-compatible
    gradient function for the duration of the test session.
    """
    with patch("numpy.gradient", new=matlab_gradient_1d_vectorized):
        yield


@pytest.fixture(scope="session", name="tmpdir")
def tmpdir_fixture():
    """Fixture to create a temporary directory for file output."""
    tmpdir_path = mkdtemp(prefix=f"disruptionpy-{time.strftime('%y%m%d-%H%M%S')}-")
    print(f"Using temporary directory: {tmpdir_path} for file output")
    yield tmpdir_path


@pytest.fixture(scope="module")
def test_file_path_f(request, tmpdir):
    """
    Fixture to generate file paths for test files.

    Parameters
    ----------
    request : FixtureRequest
        The request object for the current test.
    tmpdir : str
        The path to the temporary directory.

    Returns
    -------
    function
        A function that generates file paths with the specified suffix.
    """

    def inner(suffix):
        return os.path.join(tmpdir, f"{request.node.name}{suffix}")

    return inner


def skip_on_fast_execution(method):
    """Decorator to skip tests on fast execution environments."""
    if "GITHUB_ACTIONS" in os.environ:

        @pytest.mark.skip("fast execution")
        def wrapper(method):
            return method

        return wrapper
    return method
