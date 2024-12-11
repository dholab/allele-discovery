process CONVERT_AND_SORT {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path(sam)

    output:
    path("${file_label}.bam")

    script:
	file_label = file(sam).getSimpleName()
    """
    samtools view -bS ${sam} \
    | samtools sort -M -o ${file_label}.bam
    """

}

process SORT_BAM {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path bam

    output:
    path("*.sorted.bam")

    script:
	file_label = file(bam).getSimpleName()
    """
    cat ${bam} \
    | samtools sort -M -o ${file_label}.sorted.bam
    """

}

process CONVERT_TO_BAM {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path(sam)

    output:
    path("${file_label}.bam")

    script:
	file_label = file(sam).getSimpleName()
    """
    samtools view -bS ${sam} \
    -o ${file_label}.bam
    """

}

process INDEX_BAM {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path bam

    output:
    tuple path(bam), path("*.bam.bai")

    script:
    """
    samtools index ${bam}
    """

}

process FASTQ_CONVERSION {

    tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.fastq.gz")

    script:
    """
    samtools fastq ${bam} | bgzip -o ${sample_id}.fastq.gz
    """

}

process FAIDX {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path fasta

    output:
    tuple path(fasta), path("${fasta}.fai")

    script:
    """
    samtools faidx ${fasta}
    """

}

process FQIDX {

    tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path(fastq), path("${fastq}.fai")

    script:
    """
    samtools fqidx ${fastq}
    """

}
