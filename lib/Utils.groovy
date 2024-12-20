//
// This file holds a couple Groovy functions for counting amplicons and reference contigs
//

import java.nio.file.Files
import java.nio.file.Paths
import java.util.stream.Collectors

class Utils {

    public static void helpMessage(workflow, log, frontMatter) {

        log.info frontMatter
        log.info """
        Usage:

        The typical command for running `alleled` is as follows:
        nextflow run . 

        Mandatory arguments:

        Optional arguments:
        """
        .stripIndent()
    }

    public static void workflowDisplay(params, workflow, log, nextflow) {
        log.info """
                Workflow settings:
                -------------------------------------------------------------------------
                Launch directory            : ${workflow.launchDir}
                Workflow files directory    : ${workflow.projectDir}
                Run start time              : ${workflow.start}
                Run command                 : ${workflow.commandLine}
                Profile                     : ${workflow.profile}
                Nextflow version            : ${nextflow.version}
                Cleanup mode                : ${params.cleanup ?: ""}

                User-provided settings:
                -------------------------------------------------------------------------
                Minimum total reads         : ${params.min_total_reads}
                Minimum read length         : ${params.min_read_length}
                Quality reads to keep       : ${params.best_read_count}
                Max primer mismatches       : ${params.max_mismatch}

                User-provided inputs and outputs:
                -------------------------------------------------------------------------
                Input FASTQ directory       : ${params.input_dir ?: ""}
                TSV of primer pairs         : ${params.primer_tsv ?: ""}
                PBAA guide FASTA            : ${params.guide_fasta ?: ""}
                gDNA reference FASTA        : ${params.gdna_reference_fasta ?: ""}
                cDNA reference FASTA        : ${params.cdna_reference_fasta ?: ""}
                mRNA reference              : ${params.hla_mrna_reference ?: ""}
                Example CDS annotation      : ${params.hla_cds_annotation ?: ""}
                Results directory           : ${params.results}
                Email Address(es)           : ${params.email ?: ""}

                """
                .stripIndent()
    }

}
