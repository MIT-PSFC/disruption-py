import os
from contextlib import contextmanager
from typing import Iterable, Tuple


@contextmanager
def temporary_env_vars(env_var_list: Iterable[Tuple[str, str]]):
    """
    A context manager for temporarily modifying environment variables.

    This context manager allows you to temporarily change the values of multiple
    environment variables and restore their original values when exiting the context.

    Args:
        env_var_list (list of tuple): A list of tuples, where each tuple contains
            the name of an environment variable (str) and its new value (str).

    Usage:
        Example usage of this context manager:

        env_var_list = [("MY_ENV_VAR1", "new_value1"), ("MY_ENV_VAR2", "new_value2")]

        with temporary_env_vars(env_var_list):
            # Your code that uses the modified environment variables goes here
            print(os.environ["MY_ENV_VAR1"])
            print(os.environ["MY_ENV_VAR2"])

        # After the context manager exits, the environment variables will be
        # restored to their original values.

    Note:
        - Ensure that the values in the `env_var_list` are strings.

    """
    original_values = {}

    try:

        # Set the new environment variables and store their original values
        for key, value in env_var_list:
            original_values[key] = os.environ.get(key)
            os.environ[key] = value

        yield
    finally:
        # Restore the original environment variable values
        for key, value in env_var_list:
            if original_values[key] is None:
                del os.environ[key]
            else:
                os.environ[key] = original_values[key]
