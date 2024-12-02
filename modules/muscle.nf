process MAP_CLUSTERS_TO_CDNA {
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(query_id), val(query_seq), val(ref_id), val(ref_seq)

    output:
    path "${query_id}_to_${ref_id}.aln"

    script:
    query_seq_len = query_seq.toString().length()
    ref_seq_len = ref_seq.toString().length()
    """
    # Prepare the sequences in FASTA format
    sequences=">${query_id}
    ${query_seq}
    >${ref_id}
    ${ref_seq}"

    # Run MUSCLE alignment and count the number of matching nucleotides
    match_count=\$(echo "\${sequences}" \
    | muscle -maxiters 2 -quiet -clwstrict \
    | grep "^ " | grep -o "\\*" | wc -l)

    # Append the results to the output file in a tab-delimited format
    touch ${query_id}_to_${ref_id}.aln
    echo -e "${query_id}\t${ref_id}\t${query_seq_len}\t${ref_seq_len}\t\${match_count}" \
    >> "${query_id}_to_${ref_id}.aln"
    """
}

process COLLECT_BATCHES {
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path "results/???.aln"

    output:
    path "merged.aln"

    script:
    """
    touch merged.aln && \
    find results -type f -name "*.aln" -exec cat {} + >> merged.aln
    """
}

process COLLECT_MUSCLE_RESULTS {
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path "results/???.aln"

    output:
    path "merged.aln"

    script:
    """
    touch merged.aln && \
    find results -type f -name "*.aln" -exec cat {} + >> merged.aln
    """
}
