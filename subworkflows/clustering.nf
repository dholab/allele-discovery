include { RUN_PBAA } from "../modules/pacbio"
include { DEDUP_CLUSTERS; RENAME_WITH_IDS; SHARED_ANIMALS } from "../modules/bbmap"
include { MERGE_ALL_ANIMALS } from "../modules/seqkit"
include { RENAME_PUTATIVE_ALLELE_CLUSTERS } from "../modules/rename_cluster_seqs"


workflow CLUSTERING {

    take:
        ch_amplicons
        ch_indexed_guide

    main:

        RUN_PBAA (
            ch_amplicons
                .map { sample_id, fastq -> fastq }
                .combine( ch_indexed_guide ) 
        )

        DEDUP_CLUSTERS (
            RUN_PBAA.out
        )

        RENAME_WITH_IDS (
            DEDUP_CLUSTERS.out
        )

        MERGE_ALL_ANIMALS (
            RENAME_CLUSTERS
                .out
                .map { id, data -> data }
                .collect()
        )

        SHARED_ANIMALS (
            MERGE_PER_ANIMAL_CLUSTERS.out
        )

        RENAME_PUTATIVE_ALLELE_CLUSTERS (
            SHARED_ANIMALS.out
        )

    emit:
        RENAME_PUTATIVE_ALLELE_CLUSTERS.out

}
