include { GENOTYPE_AMPLICONS } from "../modules/minimap2"
include {
    CREATE_GENOTYPING_CSV;
    CREATE_GENOTYPING_PIVOT
    } from "../modules/genotype_table"

workflow GENOTYPE_REPORTING {

    take:
        ch_genotyping_fasta
        ch_amplicon_reads

    main:

        GENOTYPE_AMPLICONS (
            ch_genotyping_fasta,
            ch_amplicon_reads
        )

        CREATE_GENOTYPING_CSV (
            GENOTYPE_CCS.out.collect()
        )

        CREATE_GENOTYPING_PIVOT (
            CREATE_GENOTYPING_CSV.out
        )

}
