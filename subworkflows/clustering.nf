include { RUN_PBAA                        } from "../modules/pacbio"
include {
    DEDUP_CLUSTERS ;
    RENAME_CLUSTERS_WITH_IDS ;
    SHARED_ANIMALS
} from "../modules/bbmap"
include { MERGE_ALL_ANIMALS               } from "../modules/seqkit"
include { RENAME_PUTATIVE_ALLELE_CLUSTERS } from "../modules/rename_cluster_seqs"


workflow CLUSTERING {
    take:
    ch_amplicons
    ch_indexed_guide

    main:

    RUN_PBAA(
        ch_amplicons
            .map { barcode, fastq, index -> tuple( barcode, file(fastq), file(fastq).countFastq(), file(index) ) }
            .filter { barcode, _fastq, read_count, _index ->
                if (read_count < params.min_total_reads) {
                    println "The FASTQ for sample '$barcode' contains $read_count reads, which is less than the required minimum of ${params.min_total_reads}. '$barcode' will thus be ignored."
                }
                read_count >= params.min_total_reads
            }
            .map { barcode, fastq, _read_count, index -> tuple( barcode, file(fastq), file(index) ) }
            .combine(ch_indexed_guide)
    )

    DEDUP_CLUSTERS(
        RUN_PBAA.out
    )

    RENAME_CLUSTERS_WITH_IDS(
        DEDUP_CLUSTERS.out
    )

    MERGE_ALL_ANIMALS(
        RENAME_CLUSTERS_WITH_IDS.out.map { _id, data -> data }.collect()
    )

    SHARED_ANIMALS(
        MERGE_ALL_ANIMALS.out
    )

    RENAME_PUTATIVE_ALLELE_CLUSTERS(
        SHARED_ANIMALS.out
    )

    emit:
    per_sample_clusters = RENAME_CLUSTERS_WITH_IDS.out
    merged_clusters = RENAME_PUTATIVE_ALLELE_CLUSTERS.out
}
