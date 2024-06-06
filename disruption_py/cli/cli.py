#!/usr/bin/env python3

import argparse
import importlib

from disruption_py.cli.setup_script import setup, setup_check


def setup_lazy_load_script(
    parser: argparse.ArgumentParser,
    module_name: str,
    entry_point_function_name="main",
    get_parser_function_name="get_parser",
):

    @setup_check
    def run_command(args):
        module = importlib.import_module(module_name)
        entry_point_function = getattr(module, entry_point_function_name)
        get_parser_function = getattr(module, get_parser_function_name)
        script_parser: argparse.ArgumentParser = get_parser_function()
        args = script_parser.parse_args(args)
        entry_point_function(args)

    @setup_check
    def format_help():
        module = importlib.import_module(module_name)
        get_parser_function = getattr(module, get_parser_function_name)
        script_parser: argparse.ArgumentParser = get_parser_function()
        return script_parser.format_help()

    parser.set_defaults(func=run_command)
    parser.format_help = format_help


def main():
    parser = argparse.ArgumentParser(description="disruption_py command line interface")

    main_subparsers = parser.add_subparsers(dest="command")

    # setup command
    setup_parser = main_subparsers.add_parser(
        "setup", help="Setup additional dependencies for DisruptionPy"
    )
    setup_parser.set_defaults(func=setup)

    # run command
    run_parser = main_subparsers.add_parser(
        "run", help="Run scripts using DisruptionPy"
    )

    run_subparsers = run_parser.add_subparsers(title="run commands", dest="script")

    generate_datasets_parser = run_subparsers.add_parser(
        "generate_datasets", help="A basic interface with DisruptionPy"
    )
    setup_lazy_load_script(
        generate_datasets_parser, "disruption_py.cli.generate_datasets"
    )

    evaluate_parser = run_subparsers.add_parser(
        "evaluate", help="Evaluate accuracy of parameter methods in DisruptionPy"
    )
    setup_lazy_load_script(evaluate_parser, "disruption_py.cli.evaluate_methods")

    args, remaining_args = parser.parse_known_args()

    if args.command == "setup":
        args.func(args)
    elif args.command == "run":
        args.func(remaining_args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
