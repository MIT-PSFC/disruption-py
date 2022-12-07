import os
from treelib import Node, Tree

tree = Tree()
FOLDERS = [r"/fusion/projects/disruption_warning/software/matlab_programs",
           r"/fusion/projects/disruption_warning/software/peaking_factors_d3d/Physics-based_indicators/DIAG_parameterization",
           r"/fusion/projects/disruption_warning/peaking_factors_d3d/shared_scripts"]
ALL_SCRIPTS = []
for folder in FOLDERS:
    for file_name in os.listdir(folder):
        if file_name[-2:] == '.m':
            ALL_SCRIPTS.append(
                (os.path.join(folder, file_name), file_name[:-2]))


def search_script_for_paths(path):
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        if line[0].strip() != "%":
            for script, method in ALL_SCRIPTS:
                if method in line:
                    add_script(script, parent=method)


def add_script(path, parent=None):
    with open(path) as f:
        lines = f.readlines()
    method_info = lines[0]
    # Sometimes a line is too long, so we need to check if it is
    # If it is, we need to add the next line to the method_info
    for i in range(2, len(lines)):
        if lines[i][-4:] == '...':
            method_info += lines[i]
        else:
            break
    method_name, input_args, output_args = parse_method_info(method_info)
    tree.create_node(method_name + " | " + input_args +
                     " | " + output_args, method_name, parent=parent)
    search_script_for_paths(path)


# Return matlab function name, input arguments, and output arguments
def parse_method_info(method_info):
    i = method_info.find('[')
    j = method_info.find(']')
    output_args = [arg.strip() for arg in method_info[i+1:j].split(',')]
    i = method_info.find('(')
    j = method_info.find(')')
    input_args = [arg.strip() for arg in method_info[i+1:j].split(',')]
    i = method_info.find('=')
    j = method_info.find('(')
    method_name = method_info[i+1:j].strip()
    return method_name, input_args, output_args


if '__main__' == __name__:
    tree = Tree()
    add_script(
        r'/fusion/projects/disruption_warning/disruption-warning-db-workflow/D3D/disruption_warning_database_d3d.m')
    tree.show()
