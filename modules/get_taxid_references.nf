process get_taxid_reference_files{

    tag "${taxid}"
    publishDir "${params.results_dir}/reference_files/${taxid}/", mode: 'copy'


    input:
        val(taxid)
        val(kraken_db_library_path)
    output:
        tuple val(taxid), path("${taxid}.fa*"), optional : true

    shell:
    '''

    awk -v taxid="!{taxid}" '/^>/ {
        split($1, header, "|");
        if (header[2] == taxid) {
            flag=1;
            print $0 > output_file;
        } else {
            flag=0;
        }
    }
    flag' output_file="!{taxid}.fa" "!{kraken_db_library_path}" 

    # Check if the output file is not empty
    if [ -s "!{taxid}.fa" ]; then
        bwa index "!{taxid}.fa"
        echo "bwa index was run on !{taxid}.fa"
    else
        echo "Output file !{taxid}.fa does not exist."
    fi

    '''
}