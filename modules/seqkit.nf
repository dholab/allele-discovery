process FIND_COMPLETE_AMPLICONS {

  /* */

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  cpus 4

  input:
  tuple path(reads), val(amplicon_label), val(forward_primer), val(reverese_primer_rc)

  output:
  tuple val(sample_id), val(amplicon_label), val(forward_primer), val(reverese_primer_rc), path("${sample_id}_amplicons.fastq.gz")

  script:
  sample_id = file(reads).getSimpleName()
  """
  cat ${reads} | \
  seqkit grep \
  --threads ${task.cpus} \
  --max-mismatch ${params.max_mismatch} \
  --by-seq \
  --use-regexp \
  --pattern "{forward_primer}.*{reverse_primer_rc}"
  -o ${sample_id}_amplicons.fastq.gz
  """

}

process TRIM_ENDS_TO_PRIMERS {

  /* */

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  cpus 4

  input:
  tuple val(sample_id), val(amplicon_label), val(forward_primer), val(reverese_primer_rc), path(untrimmed_reads)

  output:
  tuple val(sample_id), path("${sample_id}*.trimmed.fastq.gz")

  script:
  """
  seqkit amplicon \
  -f -r 1:-1 \
  --forward ${forward_primer} \
  --reverse ${reverse_primer_rc} \
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

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  cpus 4

  input:
  tuple val(sample_id), val(amplicon_label), val(forward_primer), val(reverese_primer_rc), path("amplicons/*")

  output:
  path "${sample_id}.per_amplicon_stats.tsv"

  script:
  """
  seqkit stats \
  --threads ${task.cpus} \
  --all --basename --tabular \
  amplicons/*.fastq.gz > ${sample_id}.per_amplicon_stats.tsv
  """

}

process MERGE_BY_SAMPLE {

  /* */

  tag "${sample_id}"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  cpus 4

  input:
  tuple val(sample_id), path("fastqs/*")

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

process MERGE_ALL_ANIMALS {

  /* */

  tag "${sample_id}"

  errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
  maxRetries 2

  cpus 4

  input:
  path "per_animal_clusters/*"

  output:
  tuple val(sample_id), path("merged.fasta.gz")

  script:
  """
  seqkit scat \
  --find-only \
  --threads ${task.cpus} \
  per_animal_clusters/ \
  | bgzip -o merged.fasta.gz
  """

}

process MERGE_SEQS_FOR_GENOTYPING {

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 4

    input:
    path "genotyping_seqs/*"

    output:
    path "reads.fasta"

    script:
    """
    seqkit scat \
    --find-only \
    --threads ${task.cpus} \
    genotyping_seqs/ \
    -o reads.fasta
    """

}
