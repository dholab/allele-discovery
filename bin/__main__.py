#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "loguru",
# ]
# ///

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
    run.add_argument(
        "--nf_help",
        action="store_true",
        default=False,
        help="If supplied, uses Nextflow's help for the pipeline instead of running.",
    )
    run.add_argument(
        "--input_dir",
        type=str,
        required=True,
        help="Input directory.",
    )
    run.add_argument(
        "--primer_tsv",
        type=str,
        required=True,
        help="Primer TSV file.",
    )
    run.add_argument(
        "--guide_fasta",
        type=str,
        default=str(HERE.joinpath("assets/24923-vsearch-07-mhc-I.fasta")),
        help="Guide FASTA file (default: assets/24923-vsearch-07-mhc-I.fasta).",
    )
    run.add_argument(
        "--gdna_reference_fasta",
        type=str,
        required=True,
        help="gDNA reference FASTA.",
    )
    run.add_argument(
        "--cdna_reference_fasta",
        type=str,
        default=None,
        required=False,
        help="cDNA reference FASTA (default: None).",
    )
    run.add_argument(
        "--hla_mrna_reference",
        type=str,
        default=str(HERE.joinpath("assets/HLA-A-mRNA.fasta")),
        help="HLA mRNA reference FASTA (default: assets/HLA-A-mRNA.fasta).",
    )
    run.add_argument(
        "--hla_cds_annotation",
        type=str,
        default=str(HERE.joinpath("assets/HLA-A-mRNA-annotation.txt")),
        help="HLA CDS annotation (default: assets/HLA-A-mRNA-annotation.txt).",
    )
    run.add_argument(
        "--results",
        type=str,
        default=str(Path.cwd().joinpath("results")),
        help="Results directory (default: results in the launch directory).",
    )
    run.add_argument(
        "--min_total_reads",
        type=int,
        default=100,
        help="Minimum total reads (default: 100).",
    )
    run.add_argument(
        "--min_read_length",
        type=int,
        default=500,
        help="Minimum read length (default: 500).",
    )
    run.add_argument(
        "--best_read_count",
        type=int,
        default=1000,
        help="Best read count (default: 1000).",
    )
    run.add_argument(
        "--max_mismatch",
        type=int,
        default=0,
        help="Maximum mismatch (default: 0).",
    )
    run.add_argument(
        "--email",
        type=str,
        default=None,
        required=False,
        help="Email for pipeline notifications (default: None).",
    )
    run.add_argument(
        "--cleanup",
        action="store_true",
        default=False,
        help="If set, the pipeline performs cleanup after completion (default: False).",
    )
    run.add_argument(
        "--verbose",
        "-v",
        action="count",
        default=0,
        help="Increase verbosity level. Can be repeated, e.g. -v, -vv, -vvv.",
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

    if args.subcommands == "run":
        # If help flag is set, show Nextflow's pipeline help
        if args.nf_help:
            command = f"nextflow run {args.nf_main} --help"
            logger.info("Showing Nextflow help for the pipeline...")
        else:
            # Build the nextflow run command
            command = f"nextflow run {args.nf_main}"

            # bring in required params and params with defaults
            command += f" --input_dir {args.input_dir}"
            command += f" --results {args.results}"
            command += f" --primer_tsv {args.primer_tsv}"
            command += f" --guide_fasta {args.guide_fasta}"
            command += f" --gdna_reference_fasta {args.gdna_reference_fasta}"
            command += f" --hla_mrna_reference {args.hla_mrna_reference}"
            command += f" --hla_cds_annotation {args.hla_cds_annotation}"
            command += f" --min_total_reads {args.min_total_reads}"
            command += f" --min_read_length {args.min_read_length}"
            command += f" --best_read_count {args.best_read_count}"
            command += f" --max_mismatch {args.max_mismatch}"

            # bring in optional params
            if args.cdna_reference_fasta is not None:
                command += f" --cdna_reference_fasta {args.cdna_reference_fasta}"
            if args.cleanup:
                command += " --cleanup"
            if args.email is not None:
                command += f" --email {args.email}"

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
