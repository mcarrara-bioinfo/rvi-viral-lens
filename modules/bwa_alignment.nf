process bwa_alignment_and_post_processing {
    /*
    * Map reads to reference
    */

    publishDir "${params.results_dir}/${sample_id}_results/", overwrite: true, mode: "copy", pattern:"*.bam*"

    input:
        tuple val(file_id), path(fastq), val(reference_fasta)

    output:
        tuple val(file_id), path("${file_id}_refA.sorted.bam*")

    script:
        """
        set -e
        set -o pipefail

        #for each fasta file, get simple file name and run bwa mem
        bwa mem ${reference_fasta} ${fastq} > ${file_id}_refA.sam

        # convert sam to bam
        samtools view -S -b ${file_id}_refA.sam -o ${file_id}_refA.bam
        
        # sort alignment by leftmost coordinates
        samtools sort ${file_id}_refA.bam -o ${file_id}_refA.sorted.bam
        
        # generate indexes
        samtools index ${file_id}_refA.sorted.bam
        """
}

// note we may need to provide the index for the ref genome
// better than generate at run time for every single sample

//bwa-mem - align to reference 
//samtools-sort - sort alignments by the leftmost coordinates
