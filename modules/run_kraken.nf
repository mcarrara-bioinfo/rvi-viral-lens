process run_kraken {

    label "kraken"

    publishDir "${params.results_dir}/${sample_id}", mode: 'copy'

    input:
        tuple val(sample_id), path(fastqs) // tuple(sample_id, [fastq_pairs])
        val(db_path) // (absolute) path to kraken DB
        path(results_dir) // path to results_dir

    output:
        tuple val(sample_id), path("*.kraken.output"), path("*.class_seqs*"), path("*.unclass_seqs*"), path("*.report.txt")

    script:
        """
        #!/bin/bash

        set -e
        set -u

        kraken2 \
        --db ${db_path} \
        --output ${sample_id}.kraken.output \
        --paired \
        --classified-out ${sample_id}.class_seqs#.fq \
        --unclassified-out ${sample_id}.unclass_seqs#.fq \
        --report ${sample_id}.report.txt \
        ${fastqs}
        """
}

