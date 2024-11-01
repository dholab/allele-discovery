process FILTER_ALIGNMENTS {

	tag "${file_label}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path sam
    each path(refseq)

    output:
    path("${file_label}.filtered.bam")

    script:
	file_label = file(sam).getSimpleName()
    """
	samtools faidx ${refseq} && \
    filter_alignments.py \
    --input_sam ${sam} \
    --reference_fasta ${refseq} \
    --output_bam ${file_label}.filtered.bam
    """

}

