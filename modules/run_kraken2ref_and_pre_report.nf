process run_kraken2ref_and_pre_report {

    /*
    Run Kraken2Ref and write a report
    */
    tag "${meta.id}"
    label "kraken2ref"

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
    kraken2r -s ${meta.id} parse_report -i ${kraken_report} -o ./ -t ${params.min_reads_for_taxid}

    ## sort reads by reference (requires parse_report to have been run before)
    kraken2r -s ${meta.id} sort_reads -fq1 ${fq1_file} -fq2 ${fq2_file} -k ${kraken_output} -r ./${meta.id}_decomposed.json -u

    k2r_report.py -i ${meta.id}_decomposed.json -r ${kraken_report} -out_suffix _pre_report.tsv
    """
}