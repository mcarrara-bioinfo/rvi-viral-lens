process bwa_alignment_and_post_processing {
    /*
    * Map reads to reference
    */

    publishDir "${params.results_dir}/", overwrite: true, mode: "copy"

    input:
        tuple val(file_id), path(fastq), path(reference_fasta)

    script:
        bam_file="${file_id}.bam"

        """
        ls
        set -e
        set -o pipefail
        bwa mem "${reference_fasta}" "${fastq}" > ${file_id}.sam
        """
}

// note we may need to provide the index for the ref genome
// better than generate at run time for every single sample

