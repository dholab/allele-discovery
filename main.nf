#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {

    // assertions to check that the required parameters are provided correctly
    assert params.input_dir :
    """
    An input directory containing files with a `.fastq.gz` extension must be provided
    with the `--input_dir` command line argument.
    """
    assert file(params.input_dir).isDirectory() :
    """
    The input directory provided with `--input_dir` does not exist or is not a directory.
    """

    assert params.gdna_reference_fasta :
    """
    A full-length genomic DNA reference FASTA must be provided with the 
    `--gdna_reference_fasta` command line argument.
    """
    assert file(params.input_dir).isFile() :
    """
    The full-length genomic DNA reference FASTA file provided with
    `--gdna_reference_fasta` does not exist.
    """

    assert params.cdna_reference_fasta :
    """
    A full-length cDNA reference FASTA must be provided with the
    `--cdna_reference_fasta` command line argument.
    """
    assert file(params.input_dir).isFile() :
    """
    The full-length cDNA reference FASTA file provided with
    `--cdna_reference_fasta` does not exist.
    """


    // input channels
    ch_fastqs = Channel
        .fromPath("${params.input_dir}/*.fastq.gz")
        .map { fq -> tuple( file(fq).getSimpleName(), file(fq) ) }

    ch_mapping_reference_fasta = Channel
        .fromPath( params.mapping_reference_fasta )
    
    ch_guide_fasta = Channel
        .fromPath( params.guide_fasta )

    ch_gdna_ref = Channel
        .fromPath ( params.gdna_reference_fasta )

    ch_cdna_ref = Channel
        .fromPath ( params.cdna_reference_fasta )

    ch_hla_mrna_ref = Channel
        .fromPath ( params.hla_mrna_reference )

    ch_hla_cds_annotation = Channel
        .fromPath ( params.hla_cds_annotation )


    // workflow declaration
    MAP_CSS_TO_REF (
        ch_mapping_reference_fasta,
        ch_fastqs
    )

    EXTRACT_SOFT_CLIPPED_CONTIGS (
        MAP_CSS_TO_REF.out
    )

    GZIP_COMPRESS_TRIMMED_FASTQ (
        EXTRACT_SOFT_CLIPPED_CONTIGS.out
    )

    INDEX_FASTQ (
        GZIP_COMPRESS_TRIMMED_FASTQ.out
    )

    RUN_PBAA (
        ch_guide_fasta,
        INDEX_FASTQ
    )

    CLUSTER_PER_SAMPLE (
        RUN_PBAA.out
    )

    RENAME_CLUSTERS (
        CLUSTER_PER_SAMPLE.out
    )

    MERGE_PER_ANIMAL_CLUSTERS (
        RENAME_CLUSTERS
            .out
            .map { id, data -> data }
            .collect()
    )

    SHARED_ANIMALS (
        MERGE_PER_ANIMAL_CLUSTERS.out
    )

    RENAME_PUTATIVE_ALLELE_CLISTERS (
        SHARED_ANIMALS.out
    )

    MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA (
        ch_gdna_ref,
        RENAME_PUTATIVE_ALLELE_CLISTERS.out
    )

    FILTER_EXACT_GDNA_MATCHES (
        MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA.out,
        RENAME_PUTATIVE_ALLELE_CLISTERS.out,
        ch_gdna_ref
    )

    MAP_SHARED_CLUSTERS_TO_CDNA (
        FILTER_EXACT_GDNA_MATCHES.out.no_gdna_match,
        ch_cdna_ref
    )

    FIND_CDNA_GDNA_MATCHES (
        MAP_SHARED_CLUSTERS_TO_CDNA.out
    )

    RENAME_CDNA_MATCHED_FASTA (
        FILTER_EXACT_GDNA_MATCHES.out.no_gdna_match,
        FIND_CDNA_GDNA_MATCHES.out
    )

    PRELIMINARY_EXONERATE (
        ch_hla_mrna_ref,
        RENAME_CDNA_MATCHED_FASTA.out,
        ch_hla_cds_annotation
    )

    PROCESS_PRELIM_GFF (
        PRELIMINARY_EXONERATE.out
    )

    PRELIM_EXONERATE_MERGE_CDS (
        PROCESS_PRELIM_GFF.out
    )

    TRIM_ANNOTATIONS (
        PRELIM_EXONERATE_MERGE_CDS.out,
        FILTER_EXACT_GDNA_MATCHES.out.no_gdna_match
    )

    EXTRACT_NOVEL_SEQUENCES (
        FILTER_EXACT_GDNA_MATCHES.out.no_gdna_match,
        RENAME_CDNA_MATCHED_FASTA.out
    )

    NOVEL_EXONERATE (
        ch_hla_mrna_ref,
        EXTRACT_NOVEL_SEQUENCES.out,
        ch_hla_cds_annotation
    )

    NOVEL_EXONERATE_PROCESS_GFFS (
        NOVEL_EXONERATE.out
    )

    NOVEL_EXONERATE_MERGE_CDS (
        NOVEL_EXONERATE_PROCESS_GFFS.out
    )

    NOVEL_TRIM_ANNOTATIONS (
        NOVEL_EXONERATE_MERGE_CDS.out,
        EXTRACT_NOVEL_SEQUENCES.out
    )

    MERGE_READS (
        ch_gdna_ref,
        EXTRACT_NOVEL_SEQUENCES.out,
        RENAME_CDNA_MATCHED_FASTA.out
    )

    CLUSTAL_ALIGN (
        MERGE_READS.out,
        EXTRACT_NOVEL_SEQUENCES.out
    )

    PARSE_DISTANCES (
        CLUSTAL_ALIGN.out.distmat,
        EXTRACT_NOVEL_SEQUENCES.out,
        RENAME_CDNA_MATCHED_FASTA.out
    )

    CREATE_GENOTYPING_FASTA (
        RENAME_PUTATIVE_ALLELE_CLISTERS.out,
        ch_gdna_ref,
        MAP_SHARED_CLUSTERS_TO_CDNA.out,
        PARSE_DISTANCES.out.novel_closest_matches
    )

    GENOTYPE_CCS (
        CREATE_GENOTYPING_FASTA.out,
        EXTRACT_SOFT_CLIPPED_CONTIGS
            .out
            .map { id, fastq -> fastq }
    )

    CREATE_GENOTYPING_CSV (
        GENOTYPE_CCS.out.collect()
    )

    CREATE_GENOTYPING_PIVOT (
        CREATE_GENOTYPING_CSV.out
    )

}


process MAP_CSS_TO_REF {

    tag "${id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 3

    input:
    each path(ref)
    tuple val(id), path(fastq)

    output:
    tuple val(id), path("${id}.sam")

    script:
    """
    minimap -t ${task.cpus} ${ref} ${fastq} -ax map-hifi --secondary=no > ${id}.sam
    """
}

process EXTRACT_SOFT_CLIPPED_CONTIGS {

    tag "${id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 1

    input:
    tuple val(id), path(sam)

    output:
    tuple val(id), path("${id}.fastq")

    script:
    """
    extract_soft_clipped_contigs.py
    """

}

process GZIP_COMPRESS_TRIMMED_FASTQ {

    tag "${id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 1

    input:
    tuple val(id), path(fastq)

    output:
    tuple val(id), path("*.fastq.gz")

    script:
    """
    gzip -c ${fastq} > ${id}.fastq.gz
    """
    
}

process INDEX_FASTQ {

    tag "${id}"
    publishDir "${params.results}/trimmed", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 1

    input:
    tuple val(id), path(fastq)

    output:
    tuple val(id), path(fastq), path("${id}.fastq.fai")

    script:
    """
    samtools faidx ${fastq}
    """

}

process RUN_PBAA {

    tag "${id}"
    publishDir "${params.results}/pbaa", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 3

    input:
    each path(guide_fasta)
    tuple val(id), path(fastq), path(fai)

    output:
    tuple val(id), path("${id}_passed_cluster_sequences.fasta")

    script:
    """
    /miniconda2/bin/pbaa cluster \
    --min-read-qv 30 \
    --max-reads-per-guide 1000 \
    --max-alignments-per-read 2000 \
    ${guide_fasta} \
    ${fastq} \
    ${id}
    """

}

process CLUSTER_PER_SAMPLE {

    tag "${id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 1

    input:
    tuple val(id), path(clustered_fasta)

    output:
    tuple val(id), path(${id}.clustered.fasta.gz)

    scripts:
    """
    dedupe.sh -Xmx1g ow=t \
    in=${clustered_fasta} \
    outbest=${id}.fasta.gz \
    fo c \
    threads=${task.cpis}
    """

}

process RENAME_CLUSTERS {

    tag "${id}"

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 1

    input:
    tuple val(id), path(clustered_fa)

    output:
    tuple val(id), path("${id}.renamed.fasta.gz")

    script:
    """
    rename.sh -Xmx1g \
    in=${clustered_fq} \
    out=${id}.renamed.fasta.gz \
    prefix=${id} \
    addprefix=t \
    threads=${task.cpus}
    """

}

process MERGE_PER_ANIMAL_CLUSTERS {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path cluster_fastas

    output:
    path "merged_clusters.fasta.gz"

    script:
    """
    zcat ${cluster_fastas} | gzip > merged_clusters.fasta.gz
    """

}

process SHARED_ANIMALS {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 1

    input:
    path(merged_clusters)

    output:
    path("putative_alleles_temp.fasta")

    scrip:
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

process RENAME_PUTATIVE_ALLELE_CLISTERS {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 1

    input:
    path(putative_alleles)

    output:
    path("putative_alleles.fasta")

    script:
    """
    #!/usr/bin/env python3

    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    # parse input FASTA
    with open("putative_alleles.fasta", "w") as handle:
        for idx, record in enumerate(SeqIO.parse(${putative_alleles}, "fasta")):
            record.id = str(current_time) + '-' + str(idx)
            SeqIO.write(record, handle, "fasta")
    """

}

process MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 3

    input:
    each path(gdna_reference_fasta)
    path(putative_alleles)

    output:
    path("all_mappings.sam")

    script:
    """
    minimap2 \
    -t ${task.cpus} \
    ${gdna_reference_fasta} \
    ${putative_alleles} \
    -ax splice -N 10000 \
    > "all_mappings.sam"
    """

}

process FILTER_EXACT_GDNA_MATCHES {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path(all_mappings)
    path(putative_alleles)
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

process MAP_SHARED_CLUSTERS_TO_CDNA {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path no_gdna_match
    path cdna_ref

    output:
    path "merged.aln"

    script:
    """
    #!/usr/bin/env python3

    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    # create gDNA sequence object
    for gdna_record in SeqIO.parse(${no_gdna_match}, "fasta"):

        # create cDNA sequence object
        for cdna_record in SeqIO.parse(${cdna_ref}, "fasta"):

            # write a tab delimited line containing the gdna_record.id, cdna_record.id, cdna_record length, and match length
            # use command substitution to get the match length
            # run MUSCLE on sequence objects, but do it in a stream to avoid a lot of file I/O
            # pipe output to CLUSTALW format
            # then count the number of '*' characters that denote matches between the two sequences
            # this works for class I
            # for class II, the cDNA can be longer than the gDNA so this doesn't work
            # if the count of natching characters equals the number of cDNA characters, write to file
            # add maxiters = 2 to accelerate processing per https://www.drive5.com/muscle/manual/compromise.html

            shell( 'echo "' + gdna_record.name  + '\t' + cdna_record.name + '\t' + str(len(gdna_record.seq)) + '\t' + str(len(cdna_record.seq)) + '\t' + ' \
                $(echo ">' + gdna_record.name + '\n' + str(gdna_record.seq) + '\n>' + cdna_record.name + '\n' + str(cdna_record.seq) + '" \
                | muscle -maxiters 2 -quiet -clwstrict 2> {output[1]} \
                | grep "^ " | grep "\*" -o | wc -l )" >> {output[0]}')
    """
}

process MAP_SHARED_CLUSTERS_TO_CDNA {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path cdna_identical

    output:
    path "matches.aln"

    shell:
    '''
    awk '{{if( $4 == $5  ) print $0}}\' !{cdna_identical} > matches.aln
    '''

}

process RENAME_CDNA_MATCHED_FASTA {

    publishDir "${params.results}/cdna-identical", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path no_gdna_match
    path matches

    output:
    path "cdna_matches.fasta"

    script:
    """
    #!/usr/bin/env python3

    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    import csv

    # read input FASTA line-by-line
    for record in SeqIO.parse(${no_gdna_match}, "fasta"):

        # parse file with gdna sequences that match cdna sequences
        with open(${matches}) as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')
            for row in reader:

                # test if name of sequence in cdna match file matches gdna sequence name
                if row[0] == record.name:
                    # update name of sequence in output file
                    record.description = row[1] + '|' + row[0]

                    # write to file
                    with open("cdna_matches.fasta", "w") as handle:
                        SeqIO.write(record, handle, "fasta")
    """

}

process PRELIMINARY_EXONERATE {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path hla_mrna_ref
    path cdna_matches
    path hla_cds_annotation

    output:
    path "mapped.gff"

    script:
    """
    exonerate \
        --showtargetgff \
        --showalignment FALSE \
        --showvulgar FALSE \
        --model cdna2genome \
        --query ${hla_mrna_ref} \
        --target ${cdna_matches} \
        --refine full \
        --annotation ${hla_cds_annotation} \
        > "mapped.gff"
    """

}

process PROCESS_PRELIM_GFF {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path gff

    output:
    path "processed.gff"

    script:
    """
    process-gff.sh \
    -e ${gff}  \
    -p $projectDir/bin/21295-exonerate_gff_to_alignment_gff3.pl \
    -o processed.gff
    """

}

process PRELIM_EXONERATE_MERGE_CDS {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path processed_gff

    output:
    path "preliminary-annotations.gff"

    shell:
    '''
    awk '{{if ($3 ~ /cds/) print $1"\t"$2"\t""CDS""\t"$4,"\t"$5"\t"$6"\t"$7"\t"$8"\t""Name=CDS;ID=CDS" }}\' !{processed_gff} >> preliminary-annotations.gff
    '''

}

process TRIM_ANNOTATIONS {

    publishDir "${params.results}/cdna-identical", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path prelim_annotations

    output:
    path "cdna.gff"

    script:
    """
    trim_annotations.py
    """

}

process EXTRACT_NOVEL_SEQUENCES {

    publishDir "${params.results}/novel", mode: 'copy', overwrite: true

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

process NOVEL_EXONERATE {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path hla_mrna_ref
    path novel_alleles
    path hla_cds_annotation

    output:
    path "mapped.gff"

    script:
    """
    exonerate \
    --showtargetgff \
    --showalignment FALSE \
    --showvulgar FALSE \
    --model cdna2genome \
    --query ${hla_mrna_ref} \
    --target ${novel_alleles} \
    --refine full \
    --annotation ${hla_cds_annotation} \
    > mapped.gff
    """
}

process NOVEL_EXONERATE_PROCESS_GFFS {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path mapped_gff

    output:
    path "processed.gff"

    script:
    """
    process-gff.sh \
    -e ${mapped_gff} \
    -p ${projectDir}/bin/21295-exonerate_gff_to_alignment_gff3.pl \
    -o processed.gff
    """

}

process NOVEL_EXONERATE_MERGE_CDS {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path processed_gff

    output:
    path "preliminary-annotations.gff"

    shell:
    '''
    awk '{{if ($3 ~ /cds/) print $1"\t"$2"\t""CDS""\t"$4,"\t"$5"\t"$6"\t"$7"\t"$8"\t""Name=CDS;ID=CDS" }}\' !{processed_gff} >> preliminary-annotations.gff
    '''

}

process NOVEL_TRIM_ANNOTATIONS {

    publishDir "${params.results}/novel", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:

    output:
    path "cdna.gff"

    script:
    """
    trim_annotations.py
    """

}

process MERGE_READS {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path gdna_ref
    path novel_alleles
    path cdna_matches

    output:
    path "reads.fasta"

    script:
    """
    cat ${gdna_ref} ${novel_alleles} ${cdna_matches} > reads.fasta
    """

}

process CLUSTAL_ALIGN {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    cpus 8

    input:
    path all_discovered
    path novel_alleles

    output:
    path "aligned.fasta", emit: clustal_aligned
    path "distances.txt", emit: distmat

    script:
    """
    touch aligned.fasta distances.txt
    clustalo \
    --infile=${all_discovered} \
    --outfile=aligned.fasta \
    --distmat-out=distances.txt \
    --threads=${task.cpus} \
    --full
    """

}

process PARSE_DISTANCES {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path distances
    path novel_alleles
    path cdna_matches

    output:
    path "novel_closest_matches.xlsx", novel_closest_matches
    path "distances_tmp.txt", distances_tmp

    script:
    """
    parse_clustalo_distances.py
    """

}

process CREATE_GENOTYPING_FASTA {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path putative_alleles
    path gdna_reference_fasta
    path matches_alignment
    path novel_closest_matches

    output:
    path "classified.fasta"

    script:
    """
    create_genotyping_fasta.py
    """

}

process GENOTYPE_CCS {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    each path(genotyping_fasta)
    path trimmed_fastq

    output:
    path "${id}.sam"

    script:
    id = file(trimmed_fastq).getSimpleName()
    """
    minimap2 \
    ${genotyping_fasta} \
    ${trimmed_fastq} \
    -ax map-hifi --eqx -t 3 \
    | reformat.sh \
    in=stdin.sam \
    out=${id}.sam \
    ref=${genotyping_fasta} \
    noheader=t \
    subfilter=0 \
    threads=1 \
    ow=t \
    -da
    """

}

process CREATE_GENOTYPING_CSV {

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path sam_files

    output:
    path "genotyping.csv"

    script:
    """
    create_genotyping_csv.py
    """

}

process CREATE_GENOTYPING_PIVOT {

    publishDir "${params.results}/genotyping", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 2

    input:
    path csv

    output:
    path "genotyping.xlsx"

    script:
    """
    genotyping.py
    """

}
