include { MERGE_SEQS_FOR_GENOTYPING } from "../modules/seqkit"
include { 
    CLUSTAL_ALIGN;
    PARSE_DISTANCES
    } from "../modules/clustalo"
include { CREATE_GENOTYPING_FASTA } from "..modules/create_genotyping_fasta"
include { GENOTYPE_AMPLICONS } from "../modules/minimap2"
include { FILTER_ALIGNMENTS } from "../modules/filter_alignments
include { SORT_BAM } from "../modules/samtools"
include { REMOVE_HEADERS } from "../modules/bbmap"

workflow PREPARE_GENOTYPES {

    take:
        ch_allele_clusters
        ch_gdna_ref
        ch_novel_seqs
        ch_cdna_matches
        ch_mapped_cdna_clusters
        ch_amplicon_reads

    main:

        MERGE_SEQS_FOR_GENOTYPING (
            ch_gdna_ref
            .mix ( ch_novel_seqs )
            .mix ( ch_cdna_matches )
            .collect()
        )

        CLUSTAL_ALIGN (
            MERGE_SEQS_FOR_GENOTYPING.out,
            ch_novel_seqs
        )

        PARSE_DISTANCES (
            CLUSTAL_ALIGN.out.distmat,
            ch_novel_seqs,
            ch_cdna_matches
        )

        CREATE_GENOTYPING_FASTA (
            ch_allele_clusters,
            ch_gdna_ref,
            ch_mapped_cdna_clusters,
            PARSE_DISTANCES.out.novel_closest_matches
        )

        GENOTYPE_AMPLICONS (
            CREATE_GENOTYPING_FASTA.out,
            ch_amplicon_reads
        )

        FILTER_ALIGNMENTS (
            GENOTYPE_AMPLICONS.out,
            CREATE_GENOTYPING_FASTA.out
        )

        SORT_BAM (
            FILTER_ALIGNMENTS.out
        )       
        
        REMOVE_HEADERS (
            SORT_BAM.out
        )

    emit:
        REMOVE_HEADERS.out
}
