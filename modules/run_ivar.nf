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


/*
samtools mpileup - The SAMtools mpileup utility provides a summary of the coverage of 
                   mapped reads on a reference sequence at a single base pair resolution

samtools mpileup: This is the main command for running the samtools mpileup tool. 
It tells the system to execute the samtools program with the mpileup subcommand.

-aa:
    This option is used to enable the generation of extended BAQ (Base Alignment Quality) values.
    BAQ is a method for estimating the likelihood that a particular base in the alignment is incorrect.
    Enabling this option may improve the accuracy of the base quality scores in the output.

-d 50000:
    This option sets the maximum read depth to 50,000.
    It limits the number of reads considered at a single genomic position. 
    This can be useful for preventing excessive memory usage when dealing 
    with very deep sequencing data, but it may also result in reduced sensitivity 
    for variant calling if the depth is too low for your specific analysis.

--reference : This option specifies the reference genome or reference sequence file that will be used for alignment and variant calling.

-a: This option is used to generate per-sample and per-base alignment statistics. It will provide additional information about the alignment, such as the number of reads covering each position and the mapping quality at each position.

-B : This option specifies the input BAM file to be processed. The BAM file should be sorted by coordinate position to work correctly with samtools.

In summary, this samtools mpileup command is used to generate a pileup of aligned sequence data from a sorted BAM file
against a reference genome or sequence. It includes options for enabling BAQ calculations, setting a maximum read depth,
generating alignment statistics, and more. 
This is a common step in bioinformatics pipelines for variant calling and other downstream analyses.
*/

// ivar

/*
ivar variants: This is the main command for running the ivar tool with the variants subcommand. 
This subcommand is typically used for variant calling in viral sequencing data.
-p : This option specifies the prefix for the output files.
-q : This option specifies the minimum mapping quality score threshold for including reads in
     the variant calling process. Reads with mapping quality scores below this threshold may be filtered out.
-t : This option sets a variant frequency threshold. In this case, it's set to 0.05, which corresponds to 5%. 
      This means that variants with a frequency below 5% in the sequencing data may be filtered out or not 
      considered as significant.
      Adjusting this threshold can impact the sensitivity and specificity of variant calling, depending on
      the specific analysis requirements.

In summary, this ivar variants command is used for variant calling on viral sequencing data. 
It allows you to specify the output file prefix, a minimum mapping quality threshold, and a
variant frequency threshold for the analysis.
*/

/*

ivar consensus: This is the main command for running the ivar tool with the consensus subcommand.
  The consensus subcommand is typically used to generate a consensus sequence from aligned
  sequencing data.
-p : This option specifies the prefix for the output files, similar to the previous command.
-q : Minimum quality score threshold to count base 
-t : This option is used to set a variant frequency threshold for including variants in the 
  consensus sequence. Setting it to 0 means that even variants with very low frequencies will be
  included in the consensus sequence.
   In other words, any variant, no matter how rare, will be considered in the consensus sequence.
-m : This option specifies the minimum depth of coverage required at a position for a base to be
  included in the consensus sequence. Positions with coverage below this threshold may have the 
  base called as 'N' (ambiguous).
-n N: This option specifies the character to be used for ambiguous or low-coverage positions in 
  the consensus sequence. In this case, it's set to 'N', which is a common choice to represent uncertainty or ambiguity in DNA sequences.

In summary, this ivar consensus command is used to generate a consensus sequence from viral sequencing data.
It allows you to specify the output file prefix, a minimum mapping quality threshold, a variant frequency
 threshold (set to 0 for inclusivity), a minimum depth threshold, and the character to use for ambiguous
 positions.

*/