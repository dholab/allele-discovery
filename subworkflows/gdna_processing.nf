include { MAP_CLUSTERS_TO_FULL_LENGTH_GDNA } from "../modules/minimap2"
include { FILTER_EXACT_GDNA_MATCHES        } from "../modules/bbmap"

workflow GDNA_PROCESSING {
    take:
    ch_gdna_ref
    ch_allele_clusters

    main:

    MAP_CLUSTERS_TO_FULL_LENGTH_GDNA(
        ch_gdna_ref,
        ch_allele_clusters
    )

    FILTER_EXACT_GDNA_MATCHES(
        MAP_CLUSTERS_TO_FULL_LENGTH_GDNA.out,
        ch_allele_clusters,
        ch_gdna_ref
    )

    emit:
    gdna_match    = FILTER_EXACT_GDNA_MATCHES.out.gdna_match
    no_gdna_match = FILTER_EXACT_GDNA_MATCHES.out.no_gdna_match
}
