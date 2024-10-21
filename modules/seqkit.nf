process FIND_COMPLETE_AMPLICONS {

  /* */

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  cpus 4

  input:
  tuple path(reads), val(amplicon_label), val(forward_primer), val(reverese_primer_rc)

  output:
  tuple val(barcode), val(amplicon_label), val(forward_primer), val(reverese_primer_rc), path("${barcode}_amplicons.fastq.gz")

  script:
  barcode = file(reads).getSimpleName()
  """
  cat ${reads} | \
  seqkit grep \
  --threads ${task.cpus} \
  --max-mismatch ${params.max_mismatch} \
  --by-seq \
  --use-regexp \
  --pattern "{forward_primer}.*{reverse_primer_rc}"
  -o ${barcode}_amplicons.fastq.gz
  """

}

process TRIM_ENDS_TO_PRIMERS {

  /* */

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  cpus 4

  input:
  tuple val(barcode), val(amplicon_label), val(forward_primer), val(reverese_primer_rc), path(untrimmed_reads)

  output:
  tuple val(barcode), path("${barcode}*.trimmed.fastq.gz")

  script:
  """
  seqkit amplicon \
  -f -r 1:-1 \
  --forward ${forward_primer} \
  --reverse ${reverse_primer_rc} \
  --max-mismatch ${params.max_mismatch} \
  --strict-mode \
  --threads ${task.cpus} \
  --out-file ${barcode}.${amplicon_label}.trimmed.fastq.gz \
  ${untrimmed_reads}
  """

}

process AMPLICON_STATS {

  /* */

  tag "${barcode}"
  publishDir params.complete_amplicons, mode: 'copy', overwrite: true

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  cpus 4

  input:
  tuple val(barcode), val(amplicon_label), val(forward_primer), val(reverese_primer_rc), path("amplicons/*")

  output:
  path "${barcode}.per_amplicon_stats.tsv"

  script:
  """
  seqkit stats \
  --threads ${task.cpus} \
  --all --basename --tabular \
  amplicons/*.fastq.gz > ${barcode}.per_amplicon_stats.tsv
  """

}

process MERGE_BY_SAMPLE {

  /* */

  tag "${barcode}"
  publishDir params.complete_amplicons, mode: 'copy', overwrite: true

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  cpus 4

  input:
  tuple val(barcode), path("fastqs/*")

  output:
  tuple val(barcode), path("${barcode}.amplicons.fastq.gz")

  script:
  """
  seqkit scat \
  --find-only \
  --threads ${task.cpus} \
  fastqs/ \
  | bgzip -o ${barcode}.amplicons.fastq.gz
  """

}

