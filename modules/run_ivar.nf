params.ivar_depth_treshold = 10 // previously 5
params.ivar_min_var_frequency = 0.25 // previously 0.60

// this process reproduce the ARTIC pipeline approach
process run_ivar{
  tag "${meta.id}"
  publishDir "${params.results_dir}/${meta.sample_id}/${meta.taxid}/", mode: "copy", pattern: "*.{fa,tsv,txt}"

  label "ivar"

  input:
    tuple val(meta), path(bams)

  output:
    tuple val(meta), path("${meta.id}.consensus.fa")

  script:
    sorted_bam = "${meta.id}.sorted.bam"
    """
    set -e
    set -o pipefail

    samtools mpileup -aa -A -B -d 0 -Q0 ${sorted_bam} | \
      ivar consensus -t ${params.ivar_min_var_frequency} -m ${params.ivar_depth_treshold} -n N -p ${meta.id}.consensus
    """
}