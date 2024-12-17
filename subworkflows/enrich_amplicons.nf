
include { DEDUPLICATE_AMPLICONS; REMOVE_SHORT_READS} from "../modules/bbmap"
include { EXTRACT_TOP_QUALITY } from "../modules/extract_top_quality"
include { ORIENT_READS } from "../modules/vsearch"
include { VALIDATE_READS } from "../modules/validate"
include {
    AMPLICON_STATS;
    FIND_COMPLETE_AMPLICONS;
    TRIM_ENDS_TO_PRIMERS;
    MERGE_BY_SAMPLE;
} from "../modules/seqkit"

workflow ENRICH_AMPLICONS {

    take:
        ch_raw_reads
        ch_primer_sets
        ch_guide_fasta

    main:

        VALIDATE_READS (
            ch_raw_reads
        )

        ORIENT_READS (
            VALIDATE_READS.out
                .map { barcode, fastq, _result -> tuple( barcode, file(fastq), file(fastq).countFastq() ) }
                .filter { barcode, _fastq, read_count ->
                    if (read_count < params.min_total_reads) {
                        println "The FASTQ for sample '$barcode' contains $read_count reads, which is less than the required minimum of ${params.min_total_reads}. '$barcode' will thus be ignored."
                    }
                    read_count >= params.min_total_reads
                }
                .map { barcode, fastq, _read_count -> tuple( barcode, file(fastq) ) },
            ch_guide_fasta
        )

        FIND_COMPLETE_AMPLICONS (
            ORIENT_READS.out
                .map { _barcode, fastq -> fastq }
                .combine ( ch_primer_sets )
        )

        TRIM_ENDS_TO_PRIMERS (
            FIND_COMPLETE_AMPLICONS.out
        )

        DEDUPLICATE_AMPLICONS (
            TRIM_ENDS_TO_PRIMERS.out
        )

        REMOVE_SHORT_READS (
            DEDUPLICATE_AMPLICONS.out
        )

        EXTRACT_TOP_QUALITY (
            REMOVE_SHORT_READS.out
        )

        AMPLICON_STATS (
            EXTRACT_TOP_QUALITY.out
        )

        MERGE_BY_SAMPLE (
            EXTRACT_TOP_QUALITY.out.groupTuple( by: 0 )
        )

    emit:
        MERGE_BY_SAMPLE.out

} 
