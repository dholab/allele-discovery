process MAP_CLUSTERS_TO_CDNA {

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple val(query_id), val(query_seq), val(ref_id), val(ref_seq)

    output:
    path "${query_id}_to_${ref_id}.aln"

    shell:
    query_seq_len = query_seq.toString().length()
    ref_seq_len = ref_seq.toString().length()
    '''
    echo "!{query_id}\t!{ref_id}\t!{query_seq_len}\t!{ref_seq_len}\t \
    $(echo '>!{query_id}\n!{query_seq}\n>!{ref_id}\n!{ref_seq}')" \
    | muscle -maxiters 2 -quiet -clwstrict \
    | grep "^ " | grep "\\*" -o | wc -l 
    '''

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
    cat results/*.aln > merged.aln
    """

}