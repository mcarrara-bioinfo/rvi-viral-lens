process run_rmRepeatFromFq {
    tag "${meta.id}"
    input:
        tuple val(meta), path(fastq_1), path(fastq_2), path(trf_out_1), path(trf_out_2)

    output:
        tuple val(meta), path("${meta.id}_trf_1.fastq"), path("${meta.id}_trf_2.fastq")

    script:
    """
    rmRepeatFromFq.py -i ${fastq_1} -t ${trf_out_1} -o ${meta.id}_trf_1.fastq
    rmRepeatFromFq.py -i ${fastq_2} -t ${trf_out_2} -o ${meta.id}_trf_2.fastq
    """
}