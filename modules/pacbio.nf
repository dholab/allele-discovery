process RUN_PBAA {
	tag "${sample_id}"
	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2
	cpus 4

	input:
	tuple val(sample_id), path(reads), path(fqidx), path(guide_fasta), path(fai)

	output:
	tuple val(sample_id), path("${sample_id}_passed_cluster_sequences.fasta")

	script:
	"""
    pbaa cluster \
    ${guide_fasta} \
    ${reads} \
    ${sample_id} \
    --skip-chimera-detection \
    --num-threads ${task.cpus}
	"""
}
