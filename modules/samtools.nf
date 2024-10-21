process CONVERT_AND_SORT {

    tag "${barcode}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(sam)

    output:
    tuple val(barcode), path("${barcode}.bam")

    script:
    """
    samtools view -bS ${sam} \
    | samtools sort -M -o ${barcode}.bam
    """

}

process SORT_BAM {

    tag "${barcode}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(bam)

    output:
    tuple val(barcode), path("${barcode}.sorted.bam")

    script:
    """
    cat ${bam} \
    | samtools sort -M -o ${barcode}.sorted.bam
    """

}

process INDEX {

    tag "${barcode}"
    publishDir params.alignment, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(bam)

    output:
    tuple val(barcode), path(bam), path("${barcode}*.bam.bai")

    script:
    """
    samtools index ${bam}
    """

}

process FASTQ_CONVERSION {

    tag "${barcode}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(bam)

    output:
    tuple val(barcode), path("${barcode}.fastq.gz")

    script:
    """
    samtools fastq ${bam} | bgzip -o ${barcode}.fastq.gz
    """

}

process FAIDX {

    tag "${barcode}"
    publishDir params.basecall_fastqs, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    tuple val(barcode), path(fastq)

    output:
    tuple val(barcode), path(fastq), path("*.fai")

    script:
    index_name = file(fastq).getName()
    """
    samtools faidx --fastq ${fastq}
    """

}
