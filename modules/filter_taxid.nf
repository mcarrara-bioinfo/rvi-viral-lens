params.min_reads_for_taxid = 100

process filter_taxids {

    input:
        tuple val(meta), path(kraken_report)

    output:
        tuple val(meta), env(taxids_lst)
    
    shell:
    '''
    #!/bin/bash
    threshold=!{params.min_reads_for_taxid}  # Replace this with your desired threshold

    taxids_lst=()  # Initialize an empty array to store taxids

    while IFS=$'\t' read -r col1 col2 reads col4 taxid col6; do
        if [[ $reads -gt $threshold ]]; then
            taxids_lst+=("$taxid")
            echo $col1 $col2 $reads $taxid $col6
        fi
    done < !{kraken_report}

    echo "${taxids_lst[@]}"  # Print the list of taxids that met the threshold
    '''
}