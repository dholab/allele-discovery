include { MAP_CLUSTERS_TO_CDNA } from ""
include { FIND_CDNA_GDNA_MATCHES } from ""
include { RENAME_CDNA_MATCHED_FASTA } from ""

workflow CDNA_PROCESSING {

    take:
        ch_no_gdna_matches
        ch_cdna_ref

    main:

        MAP_CLUSTERS_TO_CDNA (
            FILTER_EXACT_GDNA_MATCHES.out.no_gdna_match,
            ch_cdna_ref
        )

        FIND_CDNA_GDNA_MATCHES (
            MAP_SHARED_CLUSTERS_TO_CDNA.out
        )

        RENAME_CDNA_MATCHED_FASTA (
            FILTER_EXACT_GDNA_MATCHES.out.no_gdna_match,
            FIND_CDNA_GDNA_MATCHES.out
        )

        EXTRACT_NOVEL_SEQUENCES (
            ch_no_gdna_match,
            RENAME_CDNA_MATCHED_FASTA.out
        )

    emit:
        novel_seqs = EXTRACT_NOVEL_SEQUENCES.out
        no_gdna_matches = FILTER_EXACT_GDNA_MATCHES.out.no_gdna_match
        cdna_matches = RENAME_CDNA_MATCHED_FASTA.out
        mapped_cdna_clusters = MAP_CLUSTERS_TO_CDNA.out

}
