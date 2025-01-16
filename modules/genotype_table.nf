process CREATE_GENOTYPING_PIVOT {
    publishDir params.genotyping, mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path sam_files

    output:
    path "genotyping_report.xlsx"

    script:
    """
    create_genotyping_pivot.py
    """
}
