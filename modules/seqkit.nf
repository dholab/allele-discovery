process FIND_COMPLETE_AMPLICONS {

  /* */

  tag "${sample_id}, ${amplicon_label}"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  cpus 4

  input:
  tuple path(reads), val(amplicon_label), val(forward_primer), val(reverse_primer)

  output:
  tuple val(sample_id), val(amplicon_label), val(forward_primer), val(reverse_primer), path("${sample_id}.${amplicon_label}_amplicons.fastq.gz")

  script:
  sample_id = file(reads).getSimpleName()
  """
    cat ${reads} | \
    seqkit grep \
    --threads ${task.cpus} \
    --max-mismatch ${params.max_mismatch} \
    --by-seq \
    --pattern ${forward_primer},${reverse_primer} \
    -o ${sample_id}.${amplicon_label}_amplicons.fastq.gz
    """
}

process TRIM_ENDS_TO_PRIMERS {

  /* */

  tag "${sample_id}, ${amplicon_label}"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  cpus 4

  input:
  tuple val(sample_id), val(amplicon_label), val(forward_primer), val(reverse_primer), path(untrimmed_reads)

  output:
  tuple val(sample_id), path("${sample_id}*.trimmed.fastq.gz")

  script:
  fwd_len = forward_primer.length()
  rev_len = reverse_primer.length()
  """
    seqkit amplicon \
    --region ${fwd_len}:-${rev_len} \
    --forward ${forward_primer} \
    --reverse ${reverse_primer} \
    --max-mismatch ${params.max_mismatch} \
    --strict-mode \
    --threads ${task.cpus} \
    --out-file ${sample_id}.${amplicon_label}.trimmed.fastq.gz \
    ${untrimmed_reads}
    """
}

process AMPLICON_STATS {

  /* */

  tag "${sample_id}"
  publishDir params.amplicon_stats, mode: 'copy', overwrite: true

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  cpus 4

  input:
  tuple val(sample_id), path("amplicons/*")

  output:
  path "${sample_id}.per_amplicon_stats.tsv"

  script:
  """
  seqkit stats \
  --threads ${task.cpus} \
  --all --basename --tabular \
  amplicons/*.f*q* > ${sample_id}.per_amplicon_stats.tsv
  """
}

process MERGE_TOP_AMPLICONS_PER_SAMPLE {

  /* */

  tag "${sample_id}"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  cpus 4

  input:
  tuple val(sample_id), path("fastqs/????.fastq")

  output:
  tuple val(sample_id), path("${sample_id}.amplicons.fastq.gz")

  script:
  """
  seqkit scat \
  --find-only \
  --threads ${task.cpus} \
  fastqs/ \
  | bgzip -o ${sample_id}.amplicons.fastq.gz
  """
}

process MERGE_ALL_AMPLICONS_PER_SAMPLE {

  /* */

  tag "${sample_id}"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  cpus 4

  input:
  tuple val(sample_id), path("fastqs/????.fastq")

  output:
  tuple val(sample_id), path("${sample_id}.amplicons.fastq.gz")

  script:
  """
  seqkit scat \
  --find-only \
  --threads ${task.cpus} \
  fastqs/ \
  | bgzip -o ${sample_id}.amplicons.fastq.gz
  """
}

process VALIDATE_NOVEL_SEQUENCES {
  publishDir params.novel_alleles, mode: 'copy', overwrite: true

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  cpus 4

  input:
  path fasta

  output:
  path "novel.validated.fasta"

  script:
  """
    seqkit seq \
    --upper-case \
    --validate-seq \
    ${fasta} \
    -o novel.validated.fasta
    """
}

process MERGE_SEQS_FOR_GENOTYPING {

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  cpus 4

  input:
  path "genotyping_seqs/*"

  output:
  path "classification.fasta"

  script:
  """
    seqkit scat \
    --format fasta \
    --find-only \
    --threads ${task.cpus} \
    genotyping_seqs/ \
    -o classification.fasta 
    """
}
