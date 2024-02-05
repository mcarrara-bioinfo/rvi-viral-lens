process compute_depth_and_coverage {

    input:
        tuple val(meta), path(ref_fa), path(bam_file)

    output:
        tuple val(meta), env(read_depth), env(percentage_genome_coverage), env(genome_size)

    shell:
    '''
    # Perform calculations
    genome_size=$(grep -v ">" !{ref_fa} | tr -d '\n' | wc -m)
    mapped_read_count=$(samtools view -c -F 260 !{bam_file})
    read_length=$(samtools stats !{bam_file} | grep "average length" | cut -f 3)
    read_depth=$(samtools depth !{bam_file} | awk '{ sum += $3 } END { print sum/NR }')

    # if empty fasta file
    if [[ ${genome_size} -eq 0 ]]; then 
             percentage_genome_coverage=null
    else
        percentage_genome_coverage=$(((($mapped
        _read_count * $read_length) / $genome_size) * 100))
    fi
    '''
}
