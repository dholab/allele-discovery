include { RENAME_WITH_IDS } from "../modules/bbmap"
include {
    FAIDX ;
    FQIDX
} from "../modules/samtools"

workflow PREPARE_SEQUENCE_FILES {
    take:
    ch_amplicons  
    ch_guide_fasta

    main:

    RENAME_WITH_IDS(
        ch_amplicons
    )

    FQIDX(
        RENAME_WITH_IDS.out
    )

    FAIDX(
        ch_guide_fasta
    )

    emit:
    amplicons     = FQIDX.out
    indexed_guide = FAIDX.out
}
