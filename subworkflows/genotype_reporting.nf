include {
    CREATE_GENOTYPING_CSV ;
    CREATE_GENOTYPING_PIVOT ;
    NEW_PIVOT_TEST
} from "../modules/genotype_table"

workflow GENOTYPE_REPORTING {
    take:
    ch_genotyping_sams

    main:

    CREATE_GENOTYPING_CSV(
        ch_genotyping_sams.collect()
    )

    NEW_PIVOT_TEST(
        ch_genotyping_sams.collect()
    )

    CREATE_GENOTYPING_PIVOT(
        CREATE_GENOTYPING_CSV.out
    )
}
