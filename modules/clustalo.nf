process CLUSTAL_ALIGN {

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    cpus 8

    input:
    path all_discovered
    path novel_alleles

    output:
    path "aligned.fasta", emit: clustal_aligned
    path "distances.txt", emit: distmat

    script:
    """
    clustalo \
    --infile=${all_discovered} \
    --outfile=aligned.fasta \
    --distmat-out=distances.txt \
    --threads=${task.cpus} \
    --full
    """

}

process PARSE_DISTANCES {

    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
    path distances
    path novel_alleles
    path cdna_matches

    output:
    path "novel_closest_matches.xlsx", emit: novel_closest_matches
    path "distances_tmp.txt", emit: distances_tmp

    script:
    """
    parse_clustalo_distances.py \
    --novel_alleles ${novel_alleles} \
    --distances ${distances} \
    --cdna_matches ${cdna_matches}
    """

}

