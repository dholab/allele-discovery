process FILTER_ALIGNMENTS {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path bam
    each path(refseq)

    output:
    path("${file_label}.filtered.bam")

    script:
	file_label = file(bam).getSimpleName()
    """
    filter_alignments.py \
    --input_bam ${bam} \
    --reference_fasta ${refseq} \
    --output_bam ${file_label}.filtered.bam
    """

}

