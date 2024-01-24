process run_pangolin {
    publishDir "${params.results_dir}/${meta.sample_id}/${meta.taxid}/", mode: "copy", pattern: "*.csv"
    label "pangolin"
    tag "${meta.id}"    

    input:
        tuple val(meta), path(mapped_fasta)

    output:
        tuple val(meta), path(mapped_fasta), env(lineage)

    shell:
        lineage_report = "${meta.id}_lineage.csv"
        consensus_fasta = "${meta.id}.fa"
        '''
        # run pangolin
        pangolin !{consensus_fasta} --outfile !{lineage_report}
        
        # get sample lineage assignment
        # stores every value of every column as an env variable
        # PS: assumes only one row at lineage_report
        IFS=, read -r taxon lineage conflict ambiguity_score scorpio_call \
        scorpio_support scorpio_conflict scorpio_notes version \
        pangolin_version scorpio_version constellation_version \
        is_designated qc_status qc_notes note < <(tail -n +2 !{lineage_report})

        echo "Taxon: $taxon"
        echo "Lineage: $lineage"
        echo "Conflict: $conflict"
        # ... and so on for other variables
        '''
}