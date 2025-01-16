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
    echo ">${query_id}\n${query_seq}\n>${ref_id}\n${ref_seq}" \
    > sequence_pair.fasta

    # Run MUSCLE alignment and count the number of matching nucleotides
    muscle -super5 sequence_pair.fasta -output aligned_pair.fasta
    match_count=\$( cat aligned_pair.fasta \
    | seqkit range -r -1:-1 \
    | seqkit seq -s \
    | tr -d '\n-' \
    | wc -c)

    # Append the results to the output file in a tab-delimited format
    echo -e "${query_id}\t${ref_id}\t${query_seq_len}\t${ref_seq_len}\t\${match_count}" \
    > "${query_id}_to_${ref_id}.aln"
    """
}

