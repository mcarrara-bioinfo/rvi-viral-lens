process run_qc_script {
    tag {meta.id}

    publishDir "${params.results_dir}/${meta.sample_id}/${meta.taxid}/", mode: 'copy', pattern: "*.{csv,png}"

    input:
    tuple val(meta), path(bam), path(fasta), path(ref)

    output:
    tuple val(meta), path("${meta.id}.qc.csv"), path("${meta.id}.depth.png"), stdout

    script:
    qcSetting = "--illumina"
    """
    qc.py ${qcSetting}\
        --outfile ${meta.id}.qc.csv \
        --sample ${meta.id} \
        --ref ${ref} \
        --bam ${bam} \
        --fasta ${fasta} \
        --ivar_md ${params.depth_treshold}

    # print first row
    sed -n "2p" ${meta.id}.qc.csv
    """
}
