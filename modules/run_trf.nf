// those are the default option used by kneaddata
params.trf_cli_options= '2 7 7 80 10 50 500 -h -ngs' //"2 5 7 80 10 50 2000 -h -ngs"

process run_trf {
    tag "${meta.id}"
    label 'mem_1'
    label 'time_1'
    label "trf"

    input:
    tuple val(meta), path(fasta_R1), path(fasta_R2)

    output:
    tuple val(meta), path("${meta.id}_1.trf"), path("${meta.id}_2.trf"), emit: paired_trf

    script:
    """
    trf ${fasta_R1} ${params.trf_cli_options} > ${meta.id}_1.trf
    trf ${fasta_R2} ${params.trf_cli_options} > ${meta.id}_2.trf
    """
}