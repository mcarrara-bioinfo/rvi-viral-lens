process retrieve_flu_subtype_and_segment {
    tag "${meta.id}"
    input:
        tuple val(meta), val(taxid_name)

    output:
        tuple val(meta), env(flu_type), env(flu_segment) // tuple(meta, flu_type, flu_segment)

    shell:

        '''
        #!/bin/bash

        set -u
        set -o pipefail

        # Attempt to parse out the flu type and flu segment from this line
        flu_type=$(echo "!{taxid_name}" | sed -n 's/.*(\\(.*\\)).*).*\$/\\1/p')
        flu_segment=$(echo "!{taxid_name}" | grep -oP '(?<=segment ).*?(?=,)')

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
