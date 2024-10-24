include { MAP_CLUSTERS_TO_CDNA; COLLECT_MUSCLE_RESULTS } from "../modules/muscle"
include { FIND_CDNA_GDNA_MATCHES } from "../modules/awk"
include { RENAME_CDNA_MATCHED_FASTA } from "../modules/rename_cdna_matched"
include { EXTRACT_NOVEL_SEQUENCES } from "../modules/bbmap"

workflow CDNA_PROCESSING {

    take:
        ch_no_gdna_matches
        ch_cdna_ref

    main:

        MAP_CLUSTERS_TO_CDNA (
            ch_no_gdna_matches
                .splitFasta( record: [id: true, seqString: true] )
                .combine(
                    ch_cdna_ref
                        .splitFasta( record: [id: true, seqString: true] )
                )
                .map { nogdna_record, cdna_ref -> 
                    tuple( nogdna_record.id, no_gdna_record.seqString, cdna_ref.id, cdna_ref.seqString )
                }
        )

        COLLECT_MUSCLE_RESULTS (
            MAP_CLUSTERS_TO_CDNA.out.collect()
        )

        FIND_CDNA_GDNA_MATCHES (
            COLLECT_MUSCLE_RESULTS.out
        )

        RENAME_CDNA_MATCHED_FASTA (
            ch_no_gdna_matches,
            FIND_CDNA_GDNA_MATCHES.out
        )

        EXTRACT_NOVEL_SEQUENCES (
            ch_no_gdna_match,
            RENAME_CDNA_MATCHED_FASTA.out
        )

    emit:
        novel_seqs = EXTRACT_NOVEL_SEQUENCES.out
        cdna_matches = RENAME_CDNA_MATCHED_FASTA.out
        mapped_cdna_clusters = MAP_CLUSTERS_TO_CDNA.out

}
