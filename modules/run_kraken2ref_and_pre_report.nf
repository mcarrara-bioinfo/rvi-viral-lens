process run_k2r_sort_reads {
    /*
    * Parse kraken report and write tax_to_reads_json
    *
    * This process runs Kraken2Ref to parse Kraken 
    * reports and sort reads by taxonomic classification.
    * It generates JSON files that map taxonomic IDs to
    * read IDs and performs sorting based on the 
    * decomposed taxonomy tree if available.
    *
    * check docs/modules/run_kraken2ref_and_pre_report.md for more extensive documentation
    */

    tag "${meta.id} - ${task.attempt} - ${task.memory}"
    cache 'lenient'
    label 'kraken2ref'
    label 'mem_k2r_escalate'
    publishDir "${params.results_dir}/${meta.id}/reads_by_taxon/", mode: 'copy', pattern: "*.{json,tsv,log}"

    input:
        tuple val(meta), path(kraken_output), path(kraken_report)

    output:
        tuple val(meta), path("${meta.id}_tax_to_reads.json"), path("${meta.id}_decomposed.json"), optional: true, emit: json_files
    script:


    """
    kraken2ref -s ${meta.id} parse_report -i ${kraken_report} -o ./ -t ${params.min_reads_for_taxid}
    
    # if empty file, no decomposed json file will be generated
    if [ -e "${meta.id}_decomposed.json" ]; then
        kraken2ref -s ${meta.id} sort_reads -k ${kraken_output} -r ./${meta.id}_decomposed.json -m tree -u
    else
        echo "Warning: JSON file does not exist. Skipping downstream commands."
    fi

    """
}

process run_k2r_dump_fastqs_and_pre_report {
    /*
    * Write Individual Fastq Files and the Pre Report
    *
    * The run_k2r_dump_fastqs_and_pre_report process 
    * uses Kraken2Ref to extract classified reads 
    * into FASTQ files and generate a preliminary 
    * report based on the taxonomic classification 
    * data. It processes classified reads and produces
    * a detailed report for further analysis.
    *
    * check docs/modules/run_kraken2ref_and_pre_report.md for more extensive documentation
    */

    tag "${meta.id} - ${task.attempt} - ${task.memory}"
    cache 'lenient'
    label 'kraken2ref'
    memory params.k2r_dump_fq_mem//'mem_k2r_escalate'

    input:
        tuple val(meta), path(classified_fqs), path(json_tax_to_readsid), path(decomposed_json), path(kraken_report)

    output:
        tuple val(meta), path("*_R{1,2}.fq"), optional: true, emit:fq_files // tuple(meta, [id.tax_id.extracted_{1,2}.fq])
        path("${meta.id}_pre_report.tsv"), optional: true, emit: report_file //tuple (meta, log, unwritten_reads_log)

    shell:
    fq_1 = classified_fqs[0]
    fq_2 = classified_fqs[1]
    
    '''
    if [ "!{meta.splitted}" == "true"]; then
        part=$(echo !{fq_1}| awk -F'[.]' '{print $(NF-1)}')
        prefix="${part}-!{meta.id}"
    else
        prefix="!{meta.id}"
    fi

    kraken2ref -s ${prefix} dump_fastqs \
            -fq1 !{fq_1} -fq2 !{fq_2} \
            --tax_to_readsid_path !{json_tax_to_readsid} \
            -o ./ --fq_load_mode !{params.k2r_fq_load_mode} \
            -r !{decomposed_json}

    # write pre_report
    k2r_report.py -i !{decomposed_json} -r !{kraken_report} -out_suffix _pre_report.tsv
    '''
}

process concatenate_fqs_parts {
    /*
    * Concatenate splitted fastqs
    *
    * This process concatenates FASTQ files from
    * multiple parts into final combined FASTQ 
    * files for each taxonomic classification. 
    * This process ensures that all parts 
    * corresponding to the same taxonomic ID 
    * are merged into single files.
    *
    * check docs/modules/run_kraken2ref_and_pre_report.md for more extensive documentation
    */

    cache 'lenient'

    input:
        tuple val(id), path(fq_parts)

    output:
       tuple val(id), path("*_R{1,2}.fq")

    shell:

    '''
    #!/bin/bash

    # Loop through all fastq files in the current directory
    for file in *_R[12].fq; do
        # Extract taxid
        taxid=$(echo "$file"| awk -F'[_]' '{print $(NF-1)}')

        # Define the output filename
        fq_1_out="!{id}_${taxid}_R1.fq"
        fq_2_out="!{id}_${taxid}_R2.fq"

        if [ -f "$fq_1_out" ]; then
            # concatenation confirmation
            echo "Concatenation for !{id}-${taxid} already done"
        else
            # Find and concatenate all matching files into the output file
            echo "Assembing $output_filename"
            cat *-!{id}_${taxid}_R1.fq >> "$fq_1_out"
            cat *-!{id}_${taxid}_R2.fq >> "$fq_2_out"
        fi
    done
    ''' 
}