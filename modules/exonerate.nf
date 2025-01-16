process PRELIMINARY_EXONERATE {

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path hla_mrna_ref
    path cdna_matches
    path hla_cds_annotation

    output:
    path "mapped.gff"

    script:
    """
    exonerate \
        --showtargetgff \
        --showalignment FALSE \
        --showvulgar FALSE \
        --model cdna2genome \
        --query ${hla_mrna_ref} \
        --target ${cdna_matches} \
        --refine full \
        --annotation ${hla_cds_annotation} \
        > "mapped.gff"
    """

}

process NOVEL_EXONERATE {

    publishDir params.novel_annotations, mode: 'copy', overwrite: true
    
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path hla_mrna_ref
    path novel_alleles
    path hla_cds_annotation

    output:
    path "mapped.gff"

    script:
    """
    exonerate \
    --showtargetgff \
    --showalignment FALSE \
    --showvulgar FALSE \
    --model cdna2genome \
    --query ${hla_mrna_ref} \
    --target ${novel_alleles} \
    --refine full \
    --annotation ${hla_cds_annotation} \
    > mapped.gff
    """
}

