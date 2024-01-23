process generate_report_line {
    tag "${meta.id}"

    input:
        tuple val(meta), path(ref_files)

    output:
        stdout

    script:
        ref_fa = "${meta.id}.fa"
        """
        # Perform calculations
        genome_size=\$(grep -v ">" ${ref_fa} | tr -d '\n' | wc -m)
        mapped_read_count=\$(samtools view -c -F 260 ${meta.bam_file})
        read_length=\$(samtools stats ${meta.bam_file} | grep "average length" | cut -f 3)
        percentage_genome_coverage=\$((((\$mapped_read_count * \$read_length) / \$genome_size) * 100))
        read_depth=\$(samtools depth ${meta.bam_file} | awk '{ sum += \$3 } END { print sum/NR }')

        # Output final report line for this sample to stdout
        echo "${meta.id}\t${meta.virus_species}\t${meta.type}\t${meta.flu_segment}\t\$percentage_genome_coverage\%\t\$read_depth"
        """
}