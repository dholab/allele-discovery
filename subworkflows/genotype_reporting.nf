include { CREATE_GENOTYPING_PIVOT } from "../modules/genotype_table"

workflow GENOTYPE_REPORTING {
    take:
    ch_genotyping_sams

    main:

    CREATE_GENOTYPING_PIVOT(
        ch_genotyping_sams.collect()
    )
}
