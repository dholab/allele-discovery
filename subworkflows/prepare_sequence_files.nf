include { RENAME_WITH_IDS } from "../modules/bbmap"

 workflow PREPARE_SEQUENCE_FILES {

    take:
        ch_amplicons
        ch_guide_fasta

    main:

        RENAME_WITH_IDS (
            ch_amplicons
        )

        FAIDX (
            ch_guide_fasta
        )

    emit:
        amplicons = RENAME_WITH_IDS.out
        indexed_guide = FAIDX.out

}
