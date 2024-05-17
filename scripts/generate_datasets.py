from disruption_py.cli.generate_datasets import get_parser, main

if __name__ == "__main__":
    parser = get_parser()
    main(parser.parse_args())
