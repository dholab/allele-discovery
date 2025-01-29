process FILTER_ALIGNMENTS {

    tag "${file_label}"

    publishDir params.filtered_geno, mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    tuple path(sam), path(fasta), path(fai)

    output:
    path "${file_label}.filtered.sam", emit: sam
    path "*.txt", emit: stats

    script:
    file_label = file(sam).getSimpleName()
    """
    filter_alignments.py \
    --input_sam ${sam} \
    --reference_fasta ${fasta} \
    --output_sam ${file_label}.filtered.sam \
    --stats
    """
}
