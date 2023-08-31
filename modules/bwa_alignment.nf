process bwa_alignment_and_post_processing {
    /*
    * Map reads to reference
    */

    publishDir "${params.results_dir}/", overwrite: true, mode: "copy"

    input:
        tuple val(file_id), path(fastq), val(reference_fasta)

    output:
        tuple val(file_id), path("${file_id}.sam")
    script:
        """
        set -e
        set -o pipefail
        bwa mem ${reference_fasta} ${fastq} > ${file_id}.sam
        """
}

// note we may need to provide the index for the ref genome
// better than generate at run time for every single sample

