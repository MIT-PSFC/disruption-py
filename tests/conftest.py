#!/usr/bin/env python3
from tempfile import mkdtemp
from unittest.mock import patch

import pytest
import os
import shutil
import time

from disruption_py.utils.mappings.tokamak_helpers import (
    get_tokamak_from_environment,
    get_tokamak_handler,
    get_tokamak_test_columns,
    get_tokamak_test_expected_failure_columns,
    get_tokamak_test_shot_ids,
)
from disruption_py.utils.math_utils import matlab_gradient_1d_vectorized


def pytest_addoption(parser):
    parser.addoption(
        "--verbose_output", action="store_true", help="More testing information."
    )
    parser.addoption(
        "--fail_quick",
        action="store_true",
        help="Finish test and report statistics instead of failing fast.",
    )
    parser.addoption(
        "--keep-logs",
        action="store_true",
        help="Whether to delete logs directory.",
    )


@pytest.fixture(scope="session")
def verbose_output(pytestconfig):
    return pytestconfig.getoption("verbose_output")


@pytest.fixture(scope="session")
def fail_quick(pytestconfig):
    return pytestconfig.getoption("fail_quick")


@pytest.fixture(scope="session")
def keep_logs(pytestconfig):
    return pytestconfig.getoption("keep_logs")


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
@pytest.fixture(scope="session", autouse=True)
def mock_numpy_gradient():
    with patch("numpy.gradient", new=matlab_gradient_1d_vectorized):
        # The patch will be in place for the duration of the test session
        yield

@pytest.fixture(scope="session")
def tmpdir(keep_logs):
    tmpdir_path = mkdtemp(dir="tests", prefix=f"tmp-{time.strftime('%y%m%d-%H:%M:%S')}-")
    yield tmpdir_path
    if not keep_logs and os.path.exists(tmpdir_path):
        shutil.rmtree(tmpdir_path)

@pytest.fixture(scope="module")
def module_file_path_f(request, tmpdir):
    def inner(suffix):
        return os.path.join(tmpdir, f"{request.node.name}{suffix}")
    return inner

@pytest.fixture(scope="function")
def test_file_path_f(request, tmpdir):
    def inner(suffix):
        return os.path.join(tmpdir, f"{request.node.name}{suffix}")
    return inner