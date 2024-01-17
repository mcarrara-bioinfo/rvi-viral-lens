process bwa_alignment_and_post_processing {
    /*
    * Map reads to reference
    */
    tag "${meta.id}"
    publishDir "${params.results_dir}/${meta.sample_id}/${meta.taxid}/", overwrite: true, mode: "copy", pattern:"*.bam*"

    input:
        tuple val(meta), path(fastq), path(ref_files)

    output:
        tuple val(meta), path("*.bam*")

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

// note we may need to provide the index for the ref genome
// better than generate at run time for every single sample

//bwa-mem - align to reference 
//samtools-sort - sort alignments by the leftmost coordinates
