import argparse
from disruption_py.cli.generate_datasets import add_generate_datasets_arguments, generate_datasets

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate DPRF compatible datasets for training and inference. Currently only supports DIII-D data")
    add_generate_datasets_arguments(parser=parser)
    generate_datasets(parser.parse_args())
