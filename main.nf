#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREPROCESSING } from "./workflows/preprocessing"
include { ALLELE_DISCOVERY } from "./workflows/allele_discovery"
include { GENOTYPING } from "./workflows/genotyping"

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

    ch_primer_pairs = Channel
        .fromPath( params.primer_tsv )
        .splitCsv( header: true, sep: "\t", strip: true, charset: "UTF-8" )

    ch_mapping_reference_fasta = Channel
        .fromPath ( params.mapping_reference_fasta )
    
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
        PREPROCESSING.out.allele_clusters,
        PREPROCESSING.out.amplicon_reads,
        ALLELE_DISCOVERY.out.mapped_cdna_clusters,
        ALLELE_DISCOVERY.out.novel_seqs,
        ALLELE_DISCOVERY.out.cdna_matches
    )

}

