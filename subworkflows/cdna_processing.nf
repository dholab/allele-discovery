include { MAP_CLUSTERS_TO_CDNA; COLLECT_MUSCLE_RESULTS } from "../modules/muscle"
include { FIND_CDNA_GDNA_MATCHES } from "../modules/awk"
include { RENAME_CDNA_MATCHED_FASTA } from "../modules/rename_cdna_matched"
include { EXTRACT_NOVEL_SEQUENCES } from "../modules/bbmap"
include { VALIDATE_NOVEL_SEQUENCES } from "../modules/seqkit"

workflow CDNA_PROCESSING {

    take:
        ch_no_gdna_matches
        ch_cdna_ref

    main:

        if ( params.cdna_reference_fasta ) {

            ch_not_gdna_records = ch_no_gdna_matches
                .splitFasta( record: [id: true, seqString: true] )

            ch_cdna_ref_records = ch_cdna_ref
                .splitFasta( record: [id: true, seqString: true] )


            MAP_CLUSTERS_TO_CDNA (
                ch_not_gdna_records
                    .combine( ch_cdna_ref_records )
                    .map { no_gdna_record, cdna_ref -> tuple(
                        no_gdna_record.id, no_gdna_record.seqString, cdna_ref.id, cdna_ref.seqString
                        )
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
                ch_no_gdna_matches,
                RENAME_CDNA_MATCHED_FASTA.out
            )

            VALIDATE_NOVEL_SEQUENCES (
                EXTRACT_NOVEL_SEQUENCES.out
            )

        } else {
        
            VALIDATE_NOVEL_SEQUENCES (
                ch_no_gdna_matches
            )
        
        }


    emit:
        novel_seqs = VALIDATE_NOVEL_SEQUENCES.out
        cdna_matches = params.cdna_reference_fasta ? RENAME_CDNA_MATCHED_FASTA.out : Channel.empty()
        mapped_cdna_clusters = params.cdna_reference_fasta ? MAP_CLUSTERS_TO_CDNA.out : Channel.empty()

}
