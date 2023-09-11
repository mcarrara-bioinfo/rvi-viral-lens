process bwa_alignment_and_post_processing {
    /*
    * Map reads to reference
    */

    publishDir "${params.results_dir}/${file_id}_results/", overwrite: true, mode: "copy", pattern:"*.bam*"

    input:
        tuple val(file_id), val(virus_id), path(fastq), val(reference_fasta)

    output:
        tuple val(file_id), path("*.bam*")

    script:
        """
        set -e
        set -o pipefail

        #for each fasta file, get simple file name and run bwa mem
        bwa mem ${reference_fasta} ${fastq} > ${file_id}_${virus_id}.sam

        # convert sam to bam
        samtools view -S -b ${file_id}_${virus_id}.sam -o ${file_id}_${virus_id}.bam
        
        # sort alignment by leftmost coordinates
        samtools sort ${file_id}_${virus_id}.bam -o ${file_id}_${virus_id}.sorted.bam
        
        # generate indexes
        samtools index ${file_id}_${virus_id}.sorted.bam
        """
}

// note we may need to provide the index for the ref genome
// better than generate at run time for every single sample

//bwa-mem - align to reference 
//samtools-sort - sort alignments by the leftmost coordinates
