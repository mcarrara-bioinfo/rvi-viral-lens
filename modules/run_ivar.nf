params.ivar_min_depth = 10
params.ivar_freq_threshold = 0.75

// this process reproduce the ARTIC pipeline approach
process run_ivar{
  tag "${meta.id}"
  publishDir "${params.results_dir}/${meta.sample_id}/${meta.taxid}/", mode: "copy", pattern: "*.{fa,tsv,txt}"

  label "ivar"

  input:
    tuple val(meta), path(bams)

  output:
    tuple val(meta), path("${meta.id}.consensus.fa"), path(mpileup_output)

  script:
    sorted_bam = "${meta.id}.sorted.bam"
    mpileup_output="${meta.id}.mpileup.out"
    """
    set -e
    set -o pipefail

    samtools mpileup -aa -A -B -d 0 -Q0 ${sorted_bam} > ${mpileup_output}
    cat ${mpileup_output} | ivar consensus -t ${params.ivar_freq_threshold} -m ${params.ivar_min_depth} -n N -p ${meta.id}.consensus
    """
}
