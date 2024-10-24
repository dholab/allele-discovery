include {
    PRELIMINARY_EXONERATE;
    NOVEL_EXONERATE   
    } from "../modules/exonerate"
include {
    PROCESS_PRELIM_GFF;
    PRELIM_EXONERATE_MERGE_CDS;
    NOVEL_EXONERATE_PROCESS_GFF;
    NOVEL_EXONERATE_MERGE_CDS
    } from "../modules/gff_processing"
include {
    TRIM_ANNOTATIONS;
    NOVEL_TRIM_ANNOTATIONS
    } from "../modules/trim_annotations"

workflow ALLELE_ANNOTATION {

    take:
        ch_hla_mrna_ref
        ch_hla_cds_annotation
        ch_novel_seqs
        ch_no_gdna_match
        ch_cdna_matched

    main:

        PRELIMINARY_EXONERATE (
            ch_hla_mrna_ref,
            ch_cdna_matched,
            ch_hla_cds_annotation
        )

        PROCESS_PRELIM_GFF (
            PRELIMINARY_EXONERATE.out
        )

        PRELIM_EXONERATE_MERGE_CDS (
            PROCESS_PRELIM_GFF.out
        )

        TRIM_ANNOTATIONS (
            PRELIM_EXONERATE_MERGE_CDS.out,
            ch_no_gdna_match
        )

        EXTRACT_NOVEL_SEQUENCES (
            ch_no_gdna_match,
            ch_cdna_matched
        )

        NOVEL_EXONERATE (
            ch_hla_mrna_ref,
            ch_novel_seqs,
            ch_hla_cds_annotation
        )

        NOVEL_EXONERATE_PROCESS_GFFS (
            NOVEL_EXONERATE.out
        )

        NOVEL_EXONERATE_MERGE_CDS (
            NOVEL_EXONERATE_PROCESS_GFFS.out
        )

        NOVEL_TRIM_ANNOTATIONS (
            NOVEL_EXONERATE_MERGE_CDS.out,
            ch_novel_seqs
        )

}
