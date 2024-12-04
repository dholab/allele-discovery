process FIND_CDNA_GDNA_MATCHES {

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path cdna_identical

    output:
    path "matches.aln"

    script:
    """
    awk '{{if( \$4 == \$5  ) print \$0}}' ${cdna_identical} > matches.aln
    """

}
