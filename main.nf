#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { PREPROCESSING } from "./workflows/preprocessing"
include { ALLELE_DISCOVERY } from "./workflows/allele_discovery"
include { GENOTYPING } from "./workflows/genotyping"

workflow {

    frontMatter = """
        alleled: amplicon-based allele genotyping and discovery
        ================================================================================
        A bioinformatic pipeline that takes long, accurate amplicons reads and uses them
        to perform genotyping and novel allele discovery. As of 2024, the methods used
        here perform well with PacBio HiFi reads, as well as Nanopore reads generated
        with R10 chemistry and super-accuracy basecalling.
        (version 0.1.0)
        ================================================================================
        """
        .stripIndent()

    if (params.help) {
        Utils.helpMessage(workflow, log, frontMatter)
        exit 0
    }

    log.info frontMatter
    Utils.workflowDisplay(params, workflow, log, nextflow)

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

    assert params.primer_tsv :
    """
    The required input TSV of primer sequences was not provided with the --primer_tsv
    argument.
    """
    assert file(params.primer_tsv).isFile() :
    """
    The input primer TSV provided with `--primer_tsv` does not exist or is not a file.
    """

    assert params.gdna_reference_fasta :
    """
    A full-length genomic DNA reference FASTA must be provided with the 
    `--gdna_reference_fasta` command line argument.
    """
    assert file(params.gdna_reference_fasta).isFile() :
    """
    The full-length genomic DNA reference FASTA file provided with
    `--gdna_reference_fasta` does not exist.
    """

    // input channels
    ch_fastqs = Channel
        .fromPath( "${params.input_dir}/*.f*q*" )
        .map { fq -> tuple( file(fq).getSimpleName(), file(fq) ) }

    ch_primer_pairs = Channel
        .fromPath( params.primer_tsv )
        .splitCsv( header: false, sep: "\t", strip: true )
        .map { row -> tuple( row[0], row[1], row[2] ) }

    // ch_mapping_reference_fasta = Channel
    //     .fromPath ( params.mapping_reference_fasta )
    
    ch_guide_fasta = Channel
        .fromPath( params.guide_fasta )

    ch_gdna_ref = Channel
        .fromPath ( params.gdna_reference_fasta )

    ch_cdna_ref = params.cdna_reference_fasta ?
        Channel.fromPath ( params.cdna_reference_fasta ) :
        Channel.empty()

    ch_hla_mrna_ref = Channel
        .fromPath ( params.hla_mrna_reference )

    ch_hla_cds_annotation = Channel
        .fromPath ( params.hla_cds_annotation )


    PREPROCESSING (
        ch_fastqs,
        ch_primer_pairs,
        ch_guide_fasta
    )

    ALLELE_DISCOVERY (
        PREPROCESSING.out.clusters_fasta,
        ch_gdna_ref,
        ch_cdna_ref,
        ch_hla_mrna_ref,
        ch_hla_cds_annotation
    )

    GENOTYPING (
        ch_gdna_ref,
        PREPROCESSING.out.clusters_fasta,
        PREPROCESSING.out.amplicon_reads,
        ALLELE_DISCOVERY.out.mapped_cdna_clusters,
        ALLELE_DISCOVERY.out.novel_seqs,
        ALLELE_DISCOVERY.out.cdna_matches
    )

    if ( params.email ) {
        workflow.onComplete {

            def msg = """\
                allele discovery and genotyping has finished running with the following settings:

                ${Utils.workflowDisplay(params, workflow, log, nextflow)}
                """
                .stripIndent()

            sendMail(to: params.email, subject: 'alleled Execution Report', body: msg)
        }
    }

}
