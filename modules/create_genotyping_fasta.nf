process CREATE_GENOTYPING_FASTA {

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path putative_alleles
    path gdna_reference_fasta
    path matches_alignment
    path novel_closest_matches

    output:
    path "classified.fasta"

    script:
    """
    create_genotyping_fasta.py
    """

}