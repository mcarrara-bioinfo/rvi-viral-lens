process run_sra_human_scrubber {
    tag "${meta.id}"
    label "rsa_human_scrubber"
    publishDir "${params.results_dir}/${meta.id}/preprocessing/", mode: "copy"

    input:
        tuple val(meta), path(fastq_1), path(fastq_2)

    output:
        tuple val(meta), path("${meta.id}_1_clean.fastq"), path("${meta.id}_2_clean.fastq")

    script:
    // TODO check scrubber parameters to add
    """
    scrub.sh -i ${fastq_1} -o ${meta.id}_1_clean.fastq
    scrub.sh -i ${fastq_2} -o ${meta.id}_2_clean.fastq
    """
}