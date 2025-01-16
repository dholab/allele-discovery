set unstable := true

alias deps := test-deps
alias dependencies := test-deps
alias py-deps := test-py-deps

@default:
	just --list

test-sh-deps:
	@samtools version
	@exonerate --version
	@bbmap.sh -h 2> /dev/null
	@clustalo --version
	@minimap2 --version
	@gffread --version
	@pigz --version
	@picard -h
	@muscle --version
	@vsearch --version
	@perl --version
	@snakemake --version
	@pbaa --version

	echo "All dependencies have been successfully installed and are available in the command line."

[script("python3")]
test-py-deps:
	import openpyxl
	import snakemake
	import pysam
	import pandas
	import Bio
	from loguru import logger

	logger.success("All Python libraries have been successfully installed.")

test-deps: test-sh-deps test-py-deps
