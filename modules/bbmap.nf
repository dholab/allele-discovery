process CLUMP_READS {

    /* */

    tag "${sample_id}"

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

process DEDUP_CLUSTERS {

    /* */

    tag "${sample_id}"
    publishDir params.deduped_clusters, mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 2

    input:
    tuple val(sample_id), path(clusters)

    output:
    tuple val(sample_id), path("${sample_id}.clusters.deduped.fasta")

    script:
    """
    dedupe.sh -Xmx8g ow=t \
    in=${clusters} \
    outbest=${sample_id}.clusters.deduped.fasta \
    fo c \
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
    tuple val(sample_id), path("${sample_id}.amplicons.no_short.fastq")

    script:
    """
    reformat.sh \
    in=${amplicons} \
    out=${sample_id}.amplicons.no_short.fastq \
    minlength=${params.min_read_length} \
    threads=${task.cpus}
    """
}

process RENAME_WITH_IDS {

    /* */

    tag "${sample_id}"
    publishDir params.top_quality, mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 4

    input:
    tuple val(sample_id), path(sequences)

    output:
    tuple val(sample_id), path("${sample_id}.amplicons.labeled.fast*")

    script:
    output_ext = file(sequences).contains(".fasta") || file(sequences).contains(".fa")
        ? "fasta"
        : "fastq"
    """
    rename.sh -Xmx1g \
    in=${sequences} \
    out="${sample_id}.amplicons.labeled.${output_ext}" \
    prefix=${sample_id} \
    addprefix=t \
    threads=${task.cpus}
    """
}

process SHARED_ANIMALS {

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 3

    input:
    path merged_clusters

    output:
    path "putative_alleles_temp.fasta"

    script:
    """
    dedupe.sh -Xmx1g \
        in=${merged_clusters} \
        outbest=all.fasta \
        am=t ac=f arc=t fo c fcc nam=4 threads=${task.cpus}

    dedupe.sh -Xmx1g \
        in=${merged_clusters} \
        out=unique.fasta \
        am=t ac=f arc=t fo fcc uniqueonly=t threads=${task.cpus}

    dedupe.sh -Xmx1g \
        in=all.fasta,unique.fasta \
        out=shared.fasta \
        ac=f uniqueonly=t threads=${task.cpus}
    
    dedupe.sh -Xmx1g \
        in=shared.fasta \
        out=putative_alleles_temp.fasta \
        ac=t threads=${task.cpus}
    """
}


process FILTER_EXACT_GDNA_MATCHES {

    publishDir params.exact_gdna_matches, mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path all_mappings
    path putative_alleles
    path gdna_ref

    output:
    path "gdna_match.sam", emit: gdna_match
    path "no_gdna_match.fasta", emit: no_gdna_match

    script:
    """
    filterlines.sh \
	in=${all_mappings} \
	out=gdna_match.sam \
	names=NM:i:0 \
	substring=t \
	include=t
    
    filterbyname.sh \
	in=${putative_alleles} \
	names=gdna_match.sam \
	out=no_gdna_match.fasta
    """
}

process EXTRACT_NOVEL_SEQUENCES {

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 3

    input:
    path no_gdna_match
    path cdna_matches

    output:
    path "novel.fasta"

    script:
    """
    mapPacBio.sh in=${no_gdna_match} ref=${cdna_matches} outu=novel.fasta subfilter=0
    """
}

process REMOVE_HEADERS {

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 3

    input:
    path sam

    output:
    path "*.noheaders.bam"

    script:
    file_label = file(sam).getSimpleName()
    """
    reformat.sh in=${sam} out=${file_label}.noheaders.bam noheader=t -Xmx8g
    """
}
