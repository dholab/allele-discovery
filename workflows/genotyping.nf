include { PREPARE_GENOTYPES } "../subworkflows/prepare_genotypes"
include { GENOTYPE_REPORTING } "../subworkflows/genotype_reporting"

workflow GENOTYPING {

    take:
        ch_gdna_ref
        ch_allele_clusters
        ch_amplicon_reads
        ch_mapped_cdna_clusters
        ch_novel_seqs
        ch_cdna_matches

    main:

        PREPARE_GENOTYPES (
            ch_allele_clusters,
            ch_gdna_ref,
            ch_novel_seqs,
            ch_cdna_matches,
            ch_mapped_cdna_clusters
        )

        GENOTYPE_REPORTING (
            PREPARE_GENOTYPES.out,
            ch_amplicon_reads
        )

}
