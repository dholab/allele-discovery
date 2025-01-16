include { PREPARE_GENOTYPES  } from "../subworkflows/prepare_genotypes"
include { GENOTYPE_REPORTING } from "../subworkflows/genotype_reporting"

workflow GENOTYPING {
    take:
    ch_gdna_ref
    ch_allele_clusters
    ch_amplicon_clusters
    ch_mapped_cdna_clusters
    ch_novel_seqs
    ch_cdna_matches

    main:

    PREPARE_GENOTYPES(
        ch_allele_clusters,
        ch_gdna_ref,
        ch_novel_seqs,
        ch_cdna_matches,
        ch_mapped_cdna_clusters,
        ch_amplicon_clusters
    )

    GENOTYPE_REPORTING(
        PREPARE_GENOTYPES.out
    )
}
