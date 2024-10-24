include { GDNA_PROCESSING } from "../subworkflows/gdna_processing"
include { CDNA_PROCESSING } from "../subworkflows/cdna_processing"
include { ALLELE_ANNOTATION } from "../subworkflows/allele_annotation"

workflow ALLELE_DISCOVERY {

    take:
        ch_allele_clusters
        ch_gdna_ref
        ch_cdna_ref
        ch_hla_mrna_ref
        ch_hla_cds_annotation

    main:
        
        GDNA_PROCESSING (
            ch_gdna_ref,
            ch_allele_clusters
        )

        CDNA_PROCESSING (
            FILTER_EXACT_GDNA_MATCHES.out.no_gdna_match,
            ch_cdna_ref
        )

        ALLELE_ANNOTATION (
            ch_hla_mrna_ref,
            ch_hla_cds_annotation,
            DNA_PROCESSING.out.novel_seqs,
            CDNA_PROCESSING.out.no_gdna_matches,
            CDNA_PROCESSING.out.cdna_matches
        )

    emit:
        novel_seqs = CDNA_PROCESSING.out.novel_seqs
        no_gdna_matches = CDNA_PROCESSING.out.no_gdna_match
        cdna_matches = CDNA_PROCESSING.out.cdna_matches
        mapped_cdna_clusters = CDNA_PROCESSING.out.mapped_cdna_clusters

}
