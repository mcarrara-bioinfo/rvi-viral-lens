process run_kraken2ref {

    tag "${meta.id}"
    label "kraken2ref"

    publishDir "${params.results_dir}/${meta.id}/reads_by_taxon/", mode: 'copy'

    input:
        tuple val(meta), path(kraken_output), path(classified_fqs), path(kraken_report)

    output:
        tuple val(meta), path("*_R{1,2}.fq"), optional: true // tuple(meta, [id.tax_id.extracted_{1,2}.fq])


    script:
    fq1_file = classified_fqs[0]
    fq2_file = classified_fqs[1]

    """
    kraken2r -s ${meta.sample_id} parse_report -i ${kraken_report} -o ./ -t ${params.min_reads_for_taxid}

    ## sort reads by reference (requires parse_report to have been run before)
    kraken2r -s ${meta.sample_id} sort_reads -fq1 ${fq1_file} -fq2 ${fq2_file} -k ${kraken_output} -r ./kraken2ref.json -u
    """
}