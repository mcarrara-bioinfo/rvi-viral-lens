params.depth_treshold = 5
params.min_quality_score = 30
params.ivar_min_var_frequency = 0.60

// this process was based on the one available at ViralFlow (https://github.com/dezordi/ViralFlow/blob/vfnext/vfnext/modules/runIvar.nf)
process run_ivar{
  tag "${meta.id}"
  publishDir "${params.results_dir}/${meta.sample_id}/${meta.taxid}/", mode: "copy", pattern: "*.{fa,tsv,txt}"

  label "ivar"

  input:
    tuple val(meta), path(bams), path(ref_fa_fls)

  output:
    tuple val(meta), path("${consensus_out_name}.fa"), path("*.txt"), path("${meta.id}.tsv")

  script:
    sorted_bam = "${meta.id}.sorted.bam"
    depth = params.depth_treshold
    min_quality_score = params.min_quality_score
    min_var_frequency = "${params.ivar_min_var_frequency}"
    freq_suffix = min_var_frequency.replace(".", "")
    ref_fa="${meta.taxid}.fa"
    consensus_out_name = "${meta.id}.ivar${freq_suffix}"
    """
    which samtools
    samtools --version
    set -e
    set -o pipefail
    # IVAR STEP 1 ----------------------------------------------------------------
    samtools mpileup -aa -d 50000 --reference ${ref_fa} -a -B ${sorted_bam} > mpileup.out
    
    cat mpileup.out | ivar variants -p ${meta.id} -q ${min_quality_score} -t 0.05

    # IVAR STEP 2 ----------------------------------------------------------------
    cat mpileup.out | ivar consensus -p ${meta.id} -q ${min_quality_score} -t 0 -m ${depth} -n N

    # IVAR STEP 3 ----------------------------------------------------------------
    cat mpileup.out | ivar consensus -p ${consensus_out_name} -q ${min_quality_score} -t ${min_var_frequency} -n N -m ${depth}

    rm mpileup.out
    """
}