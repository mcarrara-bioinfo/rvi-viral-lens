process sort_reads_with_krakentools {

    label "kraken"

    publishDir "${params.outdir}/${sample_id}/reads_by_taxon/", mode: 'copy'

    input:
        tuple val(sample_id), path(kraken_output), path(classified_fqs), path(kraken_report)

    output:
        tuple val(sample_id), path("*.*.extracted_{1,2}.fq"), optional: true // tuple(sample_id, [sample_id.tax_id.extracted_{1,2}.fq])

    script:
        """
        #!/bin/bash

        set -e
        set -u
        set -o pipefail

        awk -F"\t" '\$3!=0 && \$4!="U"' ${kraken_report} | while read -r f1 f2 f3 f4 taxid name
            do

                extractFileBase=extract_\${taxid}
                echo -e "\${extractFileBase}\t\${taxid}\t\${name}"
                echo "running kraken_tools on \${taxid}"
                extract_kraken_reads.py \
                -k ${kraken_output} \
                -s1 ${classified_fqs[0]} \
                -s2 ${classified_fqs[1]} \
                -o ${sample_id}.\${taxid}.extracted_1.fq \
                -o2 ${sample_id}.\${taxid}.extracted_2.fq \
                -t \${taxid} \
                --report ${kraken_report} \
                --fastq-output

            done

        """
}