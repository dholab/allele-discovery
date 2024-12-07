include { WORKFLOW_INTROSPECTION        } from "../subworkflows/workflow_introspection"
include { ENRICH_AMPLICONS       } from "../subworkflows/enrich_amplicons"
include { PREPARE_SEQUENCE_FILES } from "../subworkflows/prepare_sequence_files"
include { CLUSTERING             } from "../subworkflows/clustering"


workflow PREPROCESSING {
    take:
    ch_raw_reads
    ch_primer_sets
    ch_guide_fasta

    main:

    WORKFLOW_INTROSPECTION()

    ENRICH_AMPLICONS(
        ch_raw_reads,
        ch_primer_sets,
        ch_guide_fasta
    )


    PREPARE_SEQUENCE_FILES(
        ENRICH_AMPLICONS.out,
        ch_guide_fasta
    )

    CLUSTERING(
        PREPARE_SEQUENCE_FILES.out.amplicons,
        PREPARE_SEQUENCE_FILES.out.indexed_guide
    )

    emit:
    amplicon_reads = ENRICH_AMPLICONS.out
    clusters_fasta = CLUSTERING.out
    indexed_guide  = PREPARE_SEQUENCE_FILES.out.indexed_guide
}
