import argparse
from disruption_py.cli.generate_datasets import add_generate_datasets_arguments, generate_datasets

def main():
    parser = argparse.ArgumentParser(description="disruption_py command line interface")

    subparsers = parser.add_subparsers(title="subcommands", dest="subcommand")

    # generate datasets command
    generate_datasets_command = subparsers.add_parser(
        "generate_datasets", 
        help="Generate DPRF compatible datasets for training and inference. Currently only supports CMod data."
    )
    add_generate_datasets_arguments(generate_datasets_command)
    generate_datasets_command.set_defaults(func=generate_datasets)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()