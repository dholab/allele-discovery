include {
    CREATE_GENOTYPING_CSV;
    CREATE_GENOTYPING_PIVOT
    } from "../modules/genotype_table"

workflow GENOTYPE_REPORTING {

    take:
        ch_genotyping_sams

    main:

        CREATE_GENOTYPING_CSV (
            ch_genotyping_sams.collect()
        )

        CREATE_GENOTYPING_PIVOT (
            CREATE_GENOTYPING_CSV.out
        )

}
