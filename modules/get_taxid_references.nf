process get_taxid_reference_files{

    tag "${meta.taxid}"
    publishDir "${params.results_dir}/reference_files/", mode: 'copy'


    input:
        tuple val(meta), val(kraken_db_library_path),

    output:
        tuple val(meta), path("${meta.taxid}.fa*")

    script:
    """

    # get fasta sequence from kraken db library
    awk -v \
        taxid="${meta.taxid}" -v \
        RS=">" -v \
        ORS="" '$0 ~ taxid {print ">" $0}' "${kraken_db_library_path}" \
        > "${meta.taxid}.fa"

    bwa index ${meta.taxid}
    """
}