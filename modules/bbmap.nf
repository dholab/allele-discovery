process CLUMP_READS {

    /* */

	tag "${sample_id}"
    publishDir params.merged, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), path(reads)

	output:
    tuple val(sample_id), path("${sample_id}.merged.fastq.gz")

	script:
	"""
	clumpify.sh in=${reads} out=${sample_id}.merged.fastq.gz t=${task.cpus} reorder
	"""

}

process DEDUPLICATE_AMPLICONS {

    /* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 2

	input:
	tuple val(sample_id), path(amplicons)

	output:
	tuple val(sample_id), path("${sample_id}.amplicons.deduped.fastq.gz")

	script:
	"""
    dedupe.sh \
    in=${amplicons} \
    out=${sample_id}.amplicons.deduped.fastq.gz \
    threads=${task.cpus}
    """

}


process REMOVE_SHORT_READS {

    /* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 2

	input:
	tuple val(sample_id), path(amplicons)

	output:
	tuple val(sample_id), path("${sample_id}.amplicons.no_short.fastq.gz")

	script:
	"""
    reformat.sh \
    in=${amplicons} \
    out=${sample_id}.amplicons.no_short.fastq.gz \
    minlength={params.min_read_length} \
    threads=${task.cpus}
    """

}
