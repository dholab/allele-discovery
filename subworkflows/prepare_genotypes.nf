include { MERGE_SEQS_FOR_GENOTYPING } from "../modules/seqkit"
include {
    CLUSTAL_ALIGN ;
    PARSE_DISTANCES
} from "../modules/clustalo"
include { CREATE_GENOTYPING_FASTA   } from "../modules/create_genotyping_fasta"
include { FAIDX ; CONVERT_TO_BAM ; INDEX_BAM } from "../modules/samtools"
include { GENOTYPE_AMPLICON_CLUSTERS        } from "../modules/minimap2"
include { FILTER_ALIGNMENTS         } from "../modules/filter_alignments"
include { REMOVE_HEADERS            } from "../modules/bbmap"

workflow PREPARE_GENOTYPES {
    take:
    ch_allele_clusters
    ch_gdna_ref
    ch_novel_seqs
    ch_cdna_matches
    ch_mapped_cdna_clusters
    ch_amplicon_clusters

    main:

    MERGE_SEQS_FOR_GENOTYPING(
        ch_gdna_ref.mix(ch_novel_seqs).mix(ch_cdna_matches).collect()
    )

    if (params.cdna_reference_fasta) {

        CLUSTAL_ALIGN(
            MERGE_SEQS_FOR_GENOTYPING.out,
            ch_novel_seqs
        )

        PARSE_DISTANCES(
            CLUSTAL_ALIGN.out.distmat,
            ch_novel_seqs,
            ch_cdna_matches
        )

        CREATE_GENOTYPING_FASTA(
            ch_allele_clusters,
            ch_gdna_ref,
            ch_mapped_cdna_clusters,
            PARSE_DISTANCES.out.novel_closest_matches
        )

        FAIDX(
            CREATE_GENOTYPING_FASTA.out
        )
    }
    else {

        FAIDX(
            MERGE_SEQS_FOR_GENOTYPING.out
        )
    }

    GENOTYPE_AMPLICON_CLUSTERS(
        ch_amplicon_clusters.combine(FAIDX.out)
    )

    FILTER_ALIGNMENTS(
        GENOTYPE_AMPLICON_CLUSTERS.out.combine(FAIDX.out)
    )

    emit:
    FILTER_ALIGNMENTS.out.sam
}
