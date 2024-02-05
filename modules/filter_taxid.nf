params.min_reads_for_taxid = 100

process filter_taxids {

    input:
        tuple val(meta), path(kraken_report)

    output:
        tuple val(meta), env(taxids_lst), env(taxids_lst_lvl), env(taxids_lst_names)

    shell:
    '''
    #!/bin/bash
    threshold=!{params.min_reads_for_taxid} 

    taxids_lst=()  # Initialize an empty array to store taxids
    taxids_lst_lvl=() # store taxid lvl (S2, C, G)
    taxids_lst_names=() # store taxid annotation and level

    while IFS=$'\t' read -r col1 col2 reads lvl taxid name; do
        if [[ $reads -gt $threshold ]]; then
            taxids_lst+=("$taxid")
            taxids_lst_lvl+=("$lvl")
            taxids_lst_names+=("|$name")
            echo $col1 $col2 $reads $lvl $taxid $name
        fi
    done < !{kraken_report}

    # Check if arrays are empty and set them to null if so
    [[ ${#taxids_lst[@]} -eq 0 ]] && taxids_lst=null
    [[ ${#taxids_lst_lvl[@]} -eq 0 ]] && taxids_lst_lvl=null
    [[ ${#taxids_lst_names[@]} -eq 0 ]] && taxids_lst_names=null

    echo "${taxids_lst[@]}"  # Print the list of taxids that met the threshold
    echo "${taxids_lst_lvl[@]}"
    echo "${taxids_lst_names[@]}"
    '''
}