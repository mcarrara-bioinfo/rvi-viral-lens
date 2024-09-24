process bwa_alignment_and_post_processing {
    /*
    * Map reads to reference
    *
    * The process is responsible for mapping sequencing reads
    * to a reference genome using BWA (Burrows-Wheeler Aligner),
    * followed by post-processing steps including conversion to
    * BAM format, sorting, and indexing. At current version of
    * this pipeline, this process generates high-quality, sorted
    * BAM files that are essential for consensus sequence 
    * generation.
    *
    * check `docs/bwa_alignment.md` for more extensive documentation
    */

    publishDir "${params.results_dir}/${meta.sample_id}/${meta.taxid}/", overwrite: true, mode: "copy", pattern:"*.bam*"

    input:
        tuple val(meta), path(fastq), path(ref_files)

    output:
        tuple val(meta), path("*.sorted.bam*")

    script:
        ref_fa = "${meta.taxid}.fa"
        """
        set -e
        set -o pipefail

        #for each fasta file, get simple file name and run bwa mem
        bwa mem ${ref_fa} ${fastq} > ${meta.id}.sam

        # convert sam to bam
        samtools view -S -b ${meta.id}.sam -o ${meta.id}.bam

        # sort alignment by leftmost coordinates
        samtools sort ${meta.id}.bam -o ${meta.id}.sorted.bam

        # generate indexes
        samtools index ${meta.id}.sorted.bam
        """
}
