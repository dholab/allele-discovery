process CREATE_GENOTYPING_CSV {

    publishDir params.genotyping, mode: 'copy', overwrite: true

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

    publishDir params.genotyping, mode: 'copy', overwrite: true

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path csv

    output:
    path "genotyping.xlsx"

    script:
    """
    genotyping.py ${csv} "genotyping.xlsx"
    """
}

process NEW_PIVOT_TEST {
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
