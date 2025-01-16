process RENAME_CDNA_MATCHED_FASTA {

    publishDir params.cdna_matches, mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path no_gdna_match
    path matches

    output:
    path "cdna_matches.fasta"

    script:
    """
    rename_cdna_matched.py ${no_gdna_match} ${matches}
    """
}
