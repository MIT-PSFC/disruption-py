#!/usr/bin/env python3

import os
import site

import pyodbc

from disruption_py.core.utils.enums import map_string_to_enum
from disruption_py.machine.tokamak import Tokamak
from disruption_py.machine.tokamak import get_tokamak_from_environment


def cmod_setup_check():
    """
    Run on any call to `disruption_py run ***script_name***` when using CMod.

    Ensures that the following setup has been completed:

    - The user has a logbook.sybase_login file in their home directory
    - MDSPlus is installed in the active python environment
    - The necessary ODBC Driver is installed (ODBC Driver 18 for SQL Server)
    """
    sybase_login_path = os.path.expanduser("~/logbook.sybase_login")
    if not os.path.exists(sybase_login_path):
        print(
            "setup not complete, no sybase_login file found. Please run `disruption_py setup`"
        )
        return False

    try:
        import MDSplus
    except ImportError:
        print(
            "setup not complete, MDSplus package not available. Please run `disruption_py setup`"
        )
        return False
    print("Installed: MDSplus", MDSplus.__version__)

    available_drivers = pyodbc.drivers()
    desired_driver = "ODBC Driver 18 for SQL Server"
    if desired_driver not in available_drivers:
        print(
            f"The driver '{desired_driver}' is NOT installed, this may cause issues connection to the SQL database. Please run `disruption_py setup`"
        )
    return True


def setup_check(func):
    def inner(*args):
        if len(args) != 1 or (
            not hasattr(args[0], "tokamak") or args[0].tokamak is None
        ):
            tokamak = get_tokamak_from_environment()
        else:
            tokamak = map_string_to_enum(args.tokamak, Tokamak, should_raise=False)

        if tokamak is Tokamak.CMOD:
            if not cmod_setup_check():
                return False

        return func(*args)

    return inner


def setup_cmod():
    """
    Interactive CLI program run by calling `disruption_py setup` and selecting CMod.

    Walks the user through the following setup steps:

    - Adding a logbook.sybase_login file to their home directory
    - Adding MDSPlus to the active python environment
    - Checking that the necessary ODBC Driver is installed (ODBC Driver 18 for SQL Server)
    """
    print(f"Running setup for cmod")

    # add sysbase_login
    sybase_login_path = os.path.expanduser("~/logbook.sybase_login")
    if os.path.exists(sybase_login_path):
        print("logbook.sybase_login already found in home directory, continuing setup")
    else:
        username = input("What is your username for accessing the CMOD SQL logbook? ")
        sybase_login_content = "\n".join(
            ["", "alcdb2", "logbook", f"{username}", "pfcworld"]
        )
        with open(sybase_login_path, "w") as file:
            file.write(sybase_login_content)
        print("logbook.sybase_login created in your home directory")

    # add MDSplus
    try:
        import MDSplus

        print(
            f"MDSplus {MDSplus.__version__} is already available in the current environment, continuing setup"
        )
        should_add_mdsplus = "n"
    except ImportError:
        should_add_mdsplus = input(
            "Are you on the mfe workstations and would you like to add MDSplus to your current environment? (y/n) "
        )
    if should_add_mdsplus.lower() == "y":
        site_packages_path = site.getsitepackages()[0]
        mdsplus_file_path = os.path.join(site_packages_path, "mdsplus.pth")
        mdsplus_file_content = "/usr/local/mdsplus/python"
        with open(mdsplus_file_path, "w") as file:
            file.write(mdsplus_file_content)
        print("MDSplus added to your current environment")

    # check for odbc drivers
    available_drivers = pyodbc.drivers()
    desired_driver = "ODBC Driver 18 for SQL Server"
    if desired_driver in available_drivers:
        print(
            f"The driver '{desired_driver}' is installed, access to the SQL server should be possible."
        )
    else:
        print(
            f"The driver '{desired_driver}' is NOT installed, please contact your system administrator to install it."
        )
        print(
            f"Please visit https://learn.microsoft.com/en-us/sql/connect/odbc/download-odbc-driver-for-sql-server?view=sql-server-ver16 \
            for more information on installation."
        )
        print("note: the mfe workstations have the driver installed")


def setup(*args, **kwargs):
    print(f"Welcome to the setup helper script for disruption_py")

    tokamak_string: str = input(
        'Which tokamak would you like to setup disruption_py for? ("cmod" for cmod, "d3d" for d3d, leave blank to auto-detect): '
    )
    while True:
        if tokamak_string == "":
            tokamak = get_tokamak_from_environment()
        else:
            tokamak = map_string_to_enum(tokamak_string, Tokamak, should_raise=False)

        if tokamak is None:
            tokamak_string = input(
                'Invalid input, please enter "cmod" for cmod, "d3d" for d3d, or leave blank to auto-detect: '
            )
        else:
            break

    if tokamak == Tokamak.CMOD:
        setup_cmod()
    else:
        print("This tokamak is not currently supported")
