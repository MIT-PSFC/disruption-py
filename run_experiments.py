import os
import toml

CONFIG_PATH = "./disruption_py/machine/cmod/config.toml"
INTEPRETER_PATH = "./.venv/bin/python"


# # Experimento 1
OPTION = "float32"
config = toml.load(CONFIG_PATH)
config["cmod"]["precision"]["default_precision_in"] = OPTION
config["cmod"]["precision"]["default_precision_out"] = OPTION
with open(CONFIG_PATH, "w") as config_file:
    toml.dump(config, config_file)

os.system(f"{INTEPRETER_PATH} get_data.py {OPTION}")

# # Experimento 2
OPTION = "float64"
config = toml.load(CONFIG_PATH)
config["cmod"]["precision"]["default_precision_in"] = OPTION
config["cmod"]["precision"]["default_precision_out"] = OPTION
# config["option"] = "True"
with open(CONFIG_PATH, "w") as config_file:
    toml.dump(config, config_file)

os.system(f"{INTEPRETER_PATH} get_data.py {OPTION}")

# # Experimento 3
OPTION = "None"
config = toml.load(CONFIG_PATH)
config["cmod"]["precision"]["default_precision_in"] = OPTION
config["cmod"]["precision"]["default_precision_out"] = OPTION
# config["option"] = "True"
with open(CONFIG_PATH, "w") as config_file:
    toml.dump(config, config_file)

os.system(f"{INTEPRETER_PATH} get_data.py {OPTION}")
