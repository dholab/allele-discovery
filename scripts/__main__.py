import argparse
import sys
from pathlib import Path

from loguru import logger

HERE = Path(__file__).parent


def parse_command_line_args() -> argparse.Namespace:
    """
    TODO
    """
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title="Subcommands", dest="subcommands")

    # create the run subcommand
    run = subparsers.add_parser("run", help="Run the full allele discovery pipeline")
    run.add_argument(
        "--snakefile",
        "-s",
        type=str,
        help="The location of the desired snakefile, defaulting the './Snakefile'.",
        default=HERE.joinpath(Path("Snakefile")),
    )
    run.add_argument(
        "--config",
        "-c",
        type=str,
        help="Location of the allele discovery configuration file in YAML format",
        default=HERE.joinpath(Path("config/config.yaml")),
    )
    run.add_argument(
        "--cores",
        "-j",
        type=int,
        help="The desired number of CPU cores to parallelize across.",
        default=3,
    )
    run.add_argument(
        "--dry-run",
        "-d",
        action="store_true",
        help="Whether to print the pipeline's steps without running them.",
    )
    run.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity level (-v for WARNING, -vv for INFO, -vvv for DEBUG)",
        required=False,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase verbosity level (-v for WARNING, -vv for INFO, -vvv for DEBUG)",
        required=False,
    )

    return parser.parse_args()


def main() -> None:
    args = parse_command_line_args()

    # set up logger based on the requested verbosity
    if args.verbose == 0:
        level = "WARNING"
    elif args.verbose == 1:
        level = "SUCCESS"
    elif args.verbose == 2:  # noqa: PLR2004
        level = "INFO"
    elif args.verbose == 3:  # noqa: PLR2004
        level = "DEBUG"
    else:
        level = "WARNING"
    logger.remove()
    logger.add(sys.stderr, colorize=True, level=level)

    command = (
        f"snakemake --snakefile {args.snakefile} --configfile {args.config} --cores {args.cores}"
        if not args.dry_run
        else f"snakemake --snakefile {args.snakefile} --configfile {args.config} --cores {args.cores} --dry-run"
    )

    logger.debug("Launching pipeline with the following snakemake command:")
    logger.debug(command)

    command_list = command.split(" ")
    logger.debug(f"Command prepared for subprocessing: {command_list}")

    logger.success("Pipeline successful. Goodbye.")


if __name__ == "__main__":
    main()
