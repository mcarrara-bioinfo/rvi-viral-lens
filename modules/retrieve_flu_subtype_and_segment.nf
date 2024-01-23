process retrieve_flu_subtype_and_segment {
    input:
        tuple val(meta), path(kraken_report_file) // tuple(meta, kraken_report_file)

    output:
        tuple val(meta), path(kraken_report_file), env(flu_type), env(flu_segment) // tuple(meta, updated_kraken_report_file, flu_type, flu_segment)

    shell:
        taxid = meta.taxid
        '''
        #!/bin/bash

        # Attempt to get a line from the kraken report file with this particular taxon id
        line=$(awk '\$5 == !{taxid}  { print \$0 }' !{kraken_report_file})

        # Attempt to parse out the flu type and flu segment from this line
        # Assumes that the parsing is flu specific
        flu_type=$(echo $line | sed -n 's/.*(\\(.*\\)).*).*\$/\\1/p')
        flu_segment=$(echo $line | grep -oP '(?<=segment ).*?(?=,)')

        # Set flu_type to 'Null' if the above commands return nothing
        if [[ $flu_type == "" ]]; then
            flu_type=Null
        fi
        # Set flu_segment to 'Null' if the above commands return nothing
        if [[ $flu_segment == "" ]]; then
            flu_segment=Null
        fi
        '''
}
