params.ivar_min_depth = 10
params.ivar_freq_threshold = 0.75

process run_ivar{
  /*
  * Obtain Consensus Sequences using ivar
  *
  * The run_ivar process in this Nextflow pipeline is
  * designed to reproduce the ARTIC pipeline approach
  * for generating consensus sequences from BAM files 
  * using the samtools and ivar tools. Its main 
  * objective is create consensus sequences based on
  * specified depth and frequency thresholds. The 
  * resulting files are organized by sample and 
  * taxonomic ID in the specified results directory.
  *
  * check docs/run_ivar.md for more extensive documentation
  */

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
      ivar consensus -t ${params.ivar_freq_threshold} -m ${params.ivar_min_depth} -n N -p ${meta.id}.consensus
    """
}