process run_kraken {
    /*
    * Assign Reads to Taxids
    *
    * The `run_kraken` process is designed to perform
    * taxonomic classification of sequencing reads
    * using Kraken2, a widely-used tool for rapid
    * classification of metagenomic sequences. This
    * process classifies reads into taxonomic
    * categories, and outputs classified, unclassified
    * reads, and a detailed report, organizing the
    * results by sample ID.
    *
    * check docs/run_kraken.md for more extensive documentation
    */

    tag "${meta.id}"
    label "kraken"

    publishDir "${params.results_dir}/${meta.id}", mode: 'copy'

    input:
        tuple val(meta), path(fastqs) // tuple(sample_id, [fastq_pairs])
        val(db_path) // (absolute) path to kraken DB

    output:
        tuple val(meta), path("*.kraken.output"), path("*.class_seqs*"), path("*.unclass_seqs*"), path("*.report.txt")

    script:
        """
        #!/bin/bash

        set -e
        set -u

        kraken2 \
        --db ${db_path} \
        --output ${meta.id}.kraken.output \
        --paired \
        --classified-out ${meta.id}.class_seqs#.fq \
        --unclassified-out ${meta.id}.unclass_seqs#.fq \
        --report ${meta.id}.report.txt \
        ${fastqs}
        """
}

