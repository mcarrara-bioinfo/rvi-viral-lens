params.qc_minimum_depth = 10

process run_qc_script {
    tag {meta.id}

    publishDir "${params.results_dir}/${meta.sample_id}/${meta.taxid}/", mode: 'copy', pattern: "${meta.id}.qc.csv"

    label "qc"

    input:
    tuple val(meta), path(bam), path(fasta), path(ref), path(samtools_mpileup)

    output:
    tuple val(meta), path("${meta.id}.qc.csv"), stdout

    script:
    samtools_flagstat="${meta.id}.flagstat.out"
    mpileup_depths="${meta.id}.depths.out"
    """
    # Generate required samtools flagstat file
    samtools flagstat ${bam} > ${samtools_flagstat}

    # Parse out the 3 key columns we need from mpileup output file
    cat ${samtools_mpileup} | cut -f 1,2,4 > ${mpileup_depths}

    # Run QC script
    qc.py \
        --outfile ${meta.id}.qc.csv \
        --sample ${meta.id} \
        --ref ${ref} \
        --bam ${bam} \
        --fasta ${fasta} \
        --depths_file ${mpileup_depths} \
        --flagstat_file ${samtools_flagstat} \
        --minimum_depth ${params.qc_minimum_depth} \
        --ivar_md ${params.ivar_min_depth}

    # Print first row of output file to stdout
    sed -n "2p" ${meta.id}.qc.csv
    """
}

