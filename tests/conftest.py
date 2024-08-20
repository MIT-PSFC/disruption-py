#!/usr/bin/env python3

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
    parser.addoption(
        "--verbose_output", action="store_true", help="More testing information."
    )


@pytest.fixture(scope="session")
def verbose_output(pytestconfig):
    return pytestconfig.getoption("verbose_output")


def pytest_generate_tests(metafunc):
    """Parametrize `data_column` and mark expected failure columns. Marked columns
    will xfail on assert False and xpass on assert True."""
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


@pytest.fixture(scope="session")
def tokamak():
    return resolve_tokamak_from_environment()


@pytest.fixture(scope="module")
def shotlist(tokamak):
    return get_tokamak_test_shotlist(tokamak)


@pytest.fixture(scope="module")
def data_columns(tokamak):
    return get_tokamak_test_columns(tokamak)


@pytest.fixture(scope="module")
def expected_failure_columns(tokamak):
    return get_tokamak_test_expected_failure_columns(tokamak)


# for testing against sql, values generated with matlab use a different gradient method that must be patched for testing
@pytest.fixture(scope="session", autouse=True)
def mock_numpy_gradient():
    with patch("numpy.gradient", new=matlab_gradient_1d_vectorized):
        # The patch will be in place for the duration of the test session
        yield


@pytest.fixture(scope="session")
def tmpdir():
    tmpdir_path = mkdtemp(prefix=f"disruptionpy-{time.strftime('%y%m%d-%H%M%S')}-")
    print(f"Using temporary directory: {tmpdir_path} for file output")
    yield tmpdir_path


@pytest.fixture(scope="module")
def module_file_path_f(request, tmpdir):
    def inner(suffix):
        return os.path.join(tmpdir, f"{request.node.name}{suffix}")

    return inner


@pytest.fixture(scope="module")
def test_file_path_f(request, tmpdir):
    def inner(suffix):
        return os.path.join(tmpdir, f"{request.node.name}{suffix}")

    return inner


def skip_on_fast_execution(method):
    if "GITHUB_ACTIONS" in os.environ:

        @pytest.mark.skip("fast execution")
        def wrapper(method):
            return method

        return wrapper
    return method
