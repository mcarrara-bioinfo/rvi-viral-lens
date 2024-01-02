process run_kraken {
    tag "${meta.id}"
    label "kraken"

    publishDir "${params.results_dir}/${meta.id}"//, mode: 'copy'

    input:
        tuple val(meta), path(fastqs) // tuple(sample_id, [fastq_pairs])
        val(db_path) // (absolute) path to kraken DB
        path(results_dir) // path to results_dir

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

