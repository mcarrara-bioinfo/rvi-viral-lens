process write_sorted_manifest {

    label "kraken"

    publishDir "${params.outdir}/${sample_id}", pattern: "*.csv", mode: 'copy', overwrite: true

    input:
        tuple val(sample_id), val(extracted_fqs_list)

    output:
        path("sorted_manifest.csv")

    script:
    """
        #!/bin/python3

        import os
        import csv
        from collections import defaultdict

        def create_lines(fq_file_list):
            line = defaultdict(list)
            for fq_file in fq_file_list:
                basename = fq_file.split("/")[-1]
                split_base = basename.split(".")
                if split_base[1] in line.keys():
                    line[split_base[1]].append(basename)
                else:
                    line[split_base[1]] =  [basename]
            return line

        EXTRACTED_FQS_DIR = os.path.join("${params.outdir}", "${sample_id}/reads_by_taxon/")

        clean_fq_list = "${extracted_fqs_list}".strip("[]").replace(" ", "").split(",")

        sorted_clean_fq_list = create_lines(clean_fq_list)

        if os.path.exists(os.path.join("${params.outdir}", "sorted_manifest.csv")):
            with open("sorted_manifest.csv", 'a+', newline='') as outcsv:
                appender = csv.writer(outcsv)
                for idx, tax_id in enumerate(sorted_clean_fq_list.keys()):
                appender.writerow([idx, "${sample_id}", EXTRACTED_FQS_DIR+sorted_clean_fq_list[tax_id][0], EXTRACTED_FQS_DIR+sorted_clean_fq_list[tax_id][1], tax_id])
        else:
            with open("sorted_manifest.csv", 'w', newline='') as outcsv:
                writer = csv.writer(outcsv)
                writer.writerow(["idx", "sample_id", "reads_1", "reads_2", "tax_id"])
                for idx, tax_id in enumerate(sorted_clean_fq_list.keys()):
                    writer.writerow([idx, "${sample_id}", EXTRACTED_FQS_DIR+sorted_clean_fq_list[tax_id][0], EXTRACTED_FQS_DIR+sorted_clean_fq_list[tax_id][1], tax_id])
    """
}