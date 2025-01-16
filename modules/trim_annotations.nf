process TRIM_ANNOTATIONS {

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path prelim_annotations
    path no_gdna_matches

    output:
    path "cdna.gff"

    script:
    """
    trim_annotations.py ${prelim_annotations} ${no_gdna_matches}
    """
}

process NOVEL_TRIM_ANNOTATIONS {

    publishDir params.novel_annotations, mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path prelim_gff
    path no_gdna_match

    output:
    path "cdna.gff"

    script:
    """
    trim_annotations.py ${prelim_gff} ${no_gdna_match}
    """
}
