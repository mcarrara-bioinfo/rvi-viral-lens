process run_kraken2ref_and_pre_report {

    /*
    Run Kraken2Ref and write a report
    */
    tag "${meta.id}"
    label "kraken2ref"
    cache 'lenient'
    publishDir "${params.results_dir}/${meta.id}/reads_by_taxon/", mode: 'copy', pattern: "*.{fq,json,tsv,log}"

    input:
        tuple val(meta), path(kraken_output), path(classified_fqs), path(kraken_report)

    output:
        tuple val(meta), path("*_R{1,2}.fq"), optional: true, emit:fq_files // tuple(meta, [id.tax_id.extracted_{1,2}.fq])
        path("${meta.id}_pre_report.tsv"), optional: true, emit: report_file //tuple (meta, log, unwritten_reads_log)

    script:
    fq1_file = classified_fqs[0]
    fq2_file = classified_fqs[1]

    """
    kraken2ref -s ${meta.id} parse_report -i ${kraken_report} -o ./ -t ${params.min_reads_for_taxid}

    # Check if the JSON file exists
    if [ -e "${meta.id}_decomposed.json" ]; then
        ## sort reads by reference (requires parse_report to have been run before)
        kraken2ref -s ${meta.id} sort_reads -fq1 ${fq1_file} -fq2 ${fq2_file} -k ${kraken_output} -r ./${meta.id}_decomposed.json -m tree -u

        # write pre report
        k2r_report.py -i ${meta.id}_decomposed.json -r ${kraken_report} -out_suffix _pre_report.tsv

    else
        echo "Warning: JSON file does not exist. Skipping downstream commands."
        # generate an empty pre repot report file
        touch ${meta.id}_pre_report.tsv
    fi
    """
}