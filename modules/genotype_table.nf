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
    genotyping.py ${csv}
    """

}
