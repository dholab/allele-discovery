process EXTRACT_TOP_QUALITY {

    /* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	input:
	tuple val(sample_id), path(amplicons)

    output:
    tuple val(sample_id), path("${sample_id}.top${params.best_read_count}.fastq")

    script:
    """
    extract_top_quality.py \
    --input ${amplicons} \
    --output ${sample_id}.top${params.best_read_count}.fastq \
    --log ${sample_id}_extraction.log \
    --top_n ${params.best_read_count}
    """

}
