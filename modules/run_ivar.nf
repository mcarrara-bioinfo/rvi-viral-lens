// this process was based on the one available at ViralFlow (https://github.com/dezordi/ViralFlow/blob/vfnext/vfnext/modules/runIvar.nf)
process run_ivar{
  publishDir "${params.outDir}/${sample_id}_results/", mode: "copy", pattern: "*.{fa,tsv}"

  input:
    tuple val(sample_id), path(bams), path(ref_fa)
    //path(ref_fa)

  output:
    tuple val(sample_id), path("*.depth*.fa"), path("*.txt"), path("${sample_id}.tsv")

  script:
    sorted_bam = "${bams[0].getSimpleName()}.sorted.bam"
    depth = 5 
    mapping_quality = 30

    """
    # IVAR STEP 1 ----------------------------------------------------------------
    samtools mpileup -aa -d 50000 --reference ${ref_fa} -a -B ${sorted_bam} | \
       ivar variants -p ${sample_id} -q ${mapping_quality} -t 0.05

    # IVAR STEP 2 ----------------------------------------------------------------
    samtools mpileup -aa -d 50000 --reference ${ref_fa} -a -B ${sorted_bam} | \
       ivar consensus -p ${sample_id} -q ${mapping_quality} -t 0 -m ${d} -n N

    # IVAR STEP 3 ----------------------------------------------------------------
    samtools mpileup -aa -d 50000 --reference ${ref_fa} -a -B ${sorted_bam} | \
       ivar consensus -p ${sample_id}.ivar060 -q ${mapping_quality} -t 0.60 -n N -m ${depth}

    # EDIT FILE NAMES
    mv ${sample_id}.fa ${sample_id}.depth${depth}.fa
    mv ${sample_id}.ivar060.fa ${sample_id}.depth${depth}.amb.fa
    sed -i -e 's/>.*/>${sample_id}/g' ${sample_id}.depth${depth}.fa
    sed -i -e 's/>.*/>${sample_id}/g' ${sample_id}.depth${depth}.amb.fa
    """
}