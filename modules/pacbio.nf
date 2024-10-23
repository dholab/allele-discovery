process RUN_PBAA {

    /* */

	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

	cpus 4

	input:
	tuple path(reads), path(guide_fasta), path(fai)

	output:
    tuple val(sample_id), path("${sample_id}_passed_cluster_sequences.fasta")

	script:
    sample_id = file(reads).getSimpleName()
	"""
    pbaa cluster \
    ${guide_fasta} \
    ${reads} \
    . \
    --skip-chimera-detection \
    --num-threads ${task.cpus}
	"""

}
