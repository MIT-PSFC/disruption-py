import os

from graphviz import Digraph
from treelib import Node, Tree

tree = Tree()
dot = Digraph(graph_attr={"rankdir": "LR"})
TOP_LEVEL_FOLDERS = [
    r"/fusion/projects/disruption_warning/disruption-warning-db-workflow/D3D/",
    r"/fusion/projects/disruption_warning/software/peaking_factors_d3d/Physics-based_indicators/DIAG_parameterization",
    r"/fusion/projects/disruption_warning/software/peaking_factors_d3d/shared_scripts",
]
ALL_SCRIPTS = []


def search_folder(path):
    for file_name in os.listdir(path):
        if os.path.isdir(os.path.join(path, file_name)):
            print(file_name)
            search_folder(os.path.join(path, file_name))
        elif file_name[-2:] == ".m":
            ALL_SCRIPTS.append((os.path.join(path, file_name), file_name[:-2]))
        else:
            pass


def search_script_for_paths(path, node_tag):
    with open(path) as f:
        lines = f.readlines()
    for line in lines[5:]:
        if line[0].strip() != "%":
            for script, method in ALL_SCRIPTS:
                if method in line and method not in tree:
                    add_script(script, parent=node_tag)


def add_script(path, parent=None):
    with open(path) as f:
        lines = f.readlines()
    start_line = -1
    for i in range(len(lines)):
        if "function" in lines[i]:
            start_line = i
            break
    if start_line == -1:
        return
    method_info = lines[start_line]
    # Sometimes a line is too long, so we need to check if it is
    # If it is, we need to add the next line to the method_info
    for i in range(start_line, len(lines)):
        if "..." in lines[i]:
            method_info += lines[i + 1]
        else:
            break
    method_info = method_info.replace("...", "")
    method_name, input_args, output_args = parse_method_info(method_info)
    if method_name not in tree:
        tree.create_node(
            method_name
            + "\nInputs:"
            + str(input_args)
            + "\nOutputs:"
            + str(output_args),
            method_name,
            parent=parent,
        )
        dot.node(
            method_name,
            method_name + " | " + str(input_args) + " | " + str(output_args),
        )
        search_script_for_paths(path, method_name)
    dot.edge(parent, method_name)


# Return matlab function name, input arguments, and output arguments
def parse_method_info(method_info):
    if "[" in method_info:
        i = method_info.find("[")
        j = method_info.find("]")
        output_args = [arg.strip() for arg in method_info[i + 1 : j].split(",")]
    else:
        i = method_info.find("function") + len("function")
        j = method_info.find("=")
        output_args = [method_info[i + 1 : j].strip()]
    i = method_info.find("(")
    j = method_info.find(")")
    input_args = [arg.strip() for arg in method_info[i + 1 : j].split(",")]
    i = method_info.find("=")
    j = method_info.find("(")
    method_name = method_info[i + 1 : j].strip()
    return method_name, input_args, output_args


if "__main__" == __name__:
    for folder in TOP_LEVEL_FOLDERS:
        search_folder(folder)
    tree.create_node("disruption_warning_database", "root")
    search_script_for_paths(
        r"/fusion/projects/disruption_warning/disruption-warning-db-workflow/D3D/disruption_warning_database_d3d.m",
        "root",
    )
    tree.show()
    print(dot.source)
    dot.render(view=True)
