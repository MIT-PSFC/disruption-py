import argparse
from disruption_py.cli.generate_datasets import add_generate_datasets_arguments, generate_datasets

from disruption_py.cli.setup_helper import run_setup_helper        

def main():
    parser = argparse.ArgumentParser(description="disruption_py command line interface")

    subparsers = parser.add_subparsers(title="subcommands", dest="subcommand")

    # install command
    installation_command = subparsers.add_parser("setup", help="Setup help for CMod")
    installation_command.set_defaults(func=run_setup_helper)

    # generate datasets command
    generate_datasets_command = subparsers.add_parser("generate_datasets", help="Description for command 2")
    add_generate_datasets_arguments(generate_datasets_command)
    generate_datasets_command.set_defaults(func=generate_datasets)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()