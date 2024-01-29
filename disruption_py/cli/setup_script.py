import site
import os
import pyodbc
from disruption_py.utils.mappings.mappings_helpers import map_string_to_enum
from disruption_py.utils.mappings.tokamak import Tokamak
from disruption_py.utils.mappings.tokamak_helpers import get_tokamak_from_environment

def cmod_setup_check():
    sybase_login_path = os.path.expanduser('~/logbook.sybase_login')
    if not os.path.exists(sybase_login_path):
        print("disruption_py_setup not complete, no sybase_login file found. Please run disruption_py_setup")
        return False
    
    try:
        import MDSplus
    except ImportError:
        print("disruption_py_setup not complete, MDSplus package not available. Please run disruption_py_setup")
        return False
    
    available_drivers = pyodbc.drivers()
    desired_driver = 'ODBC Driver 18 for SQL Server'
    if desired_driver not in available_drivers:
        print(f"The driver '{desired_driver}' is NOT installed, this may cause issues connection to the SQL database. Please run disruption_py_setup")
    return True
    
def setup_check(func):
    def inner(args):
        if not hasattr(args, "tokamak") or args.tokamak is None:
            tokamak = get_tokamak_from_environment()
        else:
            tokamak = map_string_to_enum(args.tokamak, Tokamak, should_raise=False)
            
        if tokamak is Tokamak.CMOD:
            if not cmod_setup_check():
                return False
        
        return func(args)
    return inner

def setup_cmod():
    print(f"Running setup for cmod")

    # add sysbase_login
    sybase_login_path = os.path.expanduser('~/logbook.sybase_login')
    if os.path.exists(sybase_login_path):
        print("logbook.sybase_login already found in home directory, continuing setup")
    else:
        username = input("What is your username for accessing the CMOD SQL logbook? ")
        sybase_login_content = "\n".join(["", "alcdb2", "logbook", f"{username}", "pfcworld"])
        with open(sybase_login_path, "w") as file:
            file.write(sybase_login_content)
        print("logbook.sybase_login created in your home directory")
    
    # add MDSplus
    try:
        import MDSplus
        print("MDSplus is already available in the current environment, continuing setup")
        should_add_mdsplus = "n"
    except ImportError:
        should_add_mdsplus = input("Are you on the mfe workstations and would you like to add MDSplus to your current environment? (y/n) ")
    if should_add_mdsplus.lower() == "y":
        site_packages_path = site.getsitepackages()[0]
        mdsplus_file_path = os.path.join(site_packages_path, "mdsplus.pth")
        mdsplus_file_content = "/usr/local/mdsplus/python"
        with open(mdsplus_file_path, "w") as file:
            file.write(mdsplus_file_content)
        print("MDSplus added to your current environment")
        
    # check for odbc drivers
    available_drivers = pyodbc.drivers()
    desired_driver = 'ODBC Driver 18 for SQL Server'
    if desired_driver in available_drivers:
        print(f"The driver '{desired_driver}' is installed, access to the SQL server should be possible.")
    else:
        print(f"The driver '{desired_driver}' is NOT installed, please contact your system administrator to install it.")
        print(f"Please visit https://learn.microsoft.com/en-us/sql/connect/odbc/download-odbc-driver-for-sql-server?view=sql-server-ver16 \
            for more information on installation.")
        print("note: the mfe workstations have the driver installed")
        
def run_setup(*args, **kwargs):
    print(f"Welcome to the setup helper script for disruption_py")
    
    tokamak_string: str = input('Which tokamak would you like to setup disruption_py for? ("cmod" for cmod, "d3d" for d3d, leave blank to auto-detect): ') 
    while True:
        if tokamak_string == "":
            tokamak = get_tokamak_from_environment()
        else:
            tokamak = map_string_to_enum(tokamak_string, Tokamak, should_raise=False)
        
        if tokamak is None:
            tokamak_string = input('Invalid input, please enter "cmod" for cmod, "d3d" for d3d, or leave blank to auto-detect: ')
        else:
            break   
        
    if tokamak == Tokamak.CMOD:
        setup_cmod()
    else:
        print("This tokamak is not currently supported")