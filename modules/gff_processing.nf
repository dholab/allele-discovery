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

process NOVEL_EXONERATE_PROCESS_GFF {

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

