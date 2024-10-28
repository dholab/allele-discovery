#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path

from loguru import logger

HERE = Path(__file__).parent.parent


def parse_command_line_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title="Subcommands", dest="subcommands")

    # create the run subcommand
    run = subparsers.add_parser("run", help="Run the full allele discovery pipeline")
    run.add_argument(
        "--nf_main",
        "-m",
        type=str,
        help="The location of the pipeline entrypoint, defaulting to main.nf in the project root directory.",
        default=HERE.joinpath(Path("main.nf")),
    )

    return parser.parse_args()


def main() -> None:
    # pull in provided command line args
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

    # construct the command based on whether a dry run was requested
    command = f"nextflow run {args.nf_main}"

    # log out the command if in debug mode
    logger.debug("Launching pipeline with the following nextflow command:")
    logger.debug(command)
    command_list = command.split(" ")
    logger.debug(f"Command prepared for subprocessing: {command_list}")

    # run the pipeline
    try:
        logger.info("The pipeline wrapper CLI is not yet operational.")
        # subprocess.run(command_list, check=True, text=True, capture_output=False)  # noqa: ERA001
    except subprocess.CalledProcessError as e:
        logger.error(
            f"Pipeline failed with the following error.\n\n```\n{e}\n```\n\nAborting.",
        )
    else:
        logger.success("Pipeline successful. Goodbye.")


if __name__ == "__main__":
    main()
