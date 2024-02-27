process run_qc_script {
    tag {meta.id}

    publishDir "${params.results_dir}/qc_plots", pattern: "${meta.id}.depth.png", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(fasta), path(ref)

    output:
    path "${meta.id}.qc.csv", emit: csv
    path "${meta.id}.depth.png"

    script:
    qcSetting = "--illumina"
    //def additional_consensus = fasta_amd.name != 'NO_FILE' ? "--fasta_amd ${fasta_amd} --ivar_amd ${params.ivarAlternativeMinDepth}" : ''
    """
    qc.py ${qcSetting}\
        --outfile ${meta.id}.qc.csv \
        --sample ${meta.id} \
        --ref ${ref} \
        --bam ${bam} \
        --fasta ${fasta} \
        --ivar_md ${params.ivarMinDepth}
        #$additional_consensus
    """
}
