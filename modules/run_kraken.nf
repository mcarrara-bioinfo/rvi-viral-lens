process run_kraken {

    label "kraken"

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
        tuple val(run_id), val(sample_id), path(fastqs) // tuple(run_id, sample_id, [fastq_pairs])
        val(db_path) // (absolute) path to kraken DB
        path(outdir) // path to outdir

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

