import site
import os
import pyodbc

def setup_cmod():
    print(f"Running setup for cmod")

    # add sysbase_login
    sybase_login_path = os.path.expanduser('~/logbook.sybase_login')
    if os.path.exists(sybase_login_path):
        print("logbook.sybase_login already found in home directory, continuing setup")
    else:
        username = input("What is your username for accessing the CMOD SQL logbook?")
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
        should_add_mdsplus = input("Are you on the mfe workstations and would you like to add MDSplus to your current environment? (y/n)")
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
    tokamak_number : str = input("Which tokamak would you like to setup disruption_py for? (1 for cmod, 2 for d3d):")
    if int(tokamak_number) == 1:
        setup_cmod()
    else:
        print("This tokamak is not currently supported")