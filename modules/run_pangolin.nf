process run_pangolin {
    /*
    * Subtype SCOV2 sequence
    *
    * The run_pangolin process in this Nextflow pipeline
    * is designed to classify SARS-CoV-2 sequences into
    * lineages using the Pangolin tool. Pangolin 
    * (Phylogenetic Assignment of Named Global Outbreak
    * LINeages) is a tool widely used in the analysis of
    * SARS-CoV-2 genomes to determine the likely lineage
    * of the virus based on the consensus sequence.
    * 
    * check docs/modules/run_pangolin.nf for a more extensive documentation
    */

    publishDir "${params.results_dir}/${meta.sample_id}/${meta.taxid}/", mode: "copy", pattern: "*.csv"
    label "pangolin"
    tag "${meta.id}"

    input:
        tuple val(meta), path(mapped_fasta)

    output:
        tuple val(meta), path(mapped_fasta), env(lineage)

    shell:
        lineage_report = "${meta.id}_lineage.csv"
        consensus_fasta = mapped_fasta
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
        '''
}