process write_sorted_manifest {

    label "kraken"
    publishDir "${params.outdir}", pattern: "*.csv", mode: 'copy', overwrite: true

    input:
        val(extracted_fqs_list)

    output:
        path("sorted_manifest.csv")

    script:
        """
        #!/bin/python3

        import os
        import csv
        from collections import defaultdict

        def create_lines(input_list):
            lines = defaultdict(list)
            for data_set in input_list:
                split_data = data_set.strip("[]").replace(" [", "").split(",")
                sample_id = split_data[0]
                EXTRACTED_FQS_DIR = os.path.join("${params.outdir}", f"{sample_id}/reads_by_taxon/")
                fq_list = split_data[1:]
                for fq_file in fq_list:
                    basename = fq_file.split("/")[-1]
                    split_base = basename.split(".")
                    k = f"{sample_id}, {split_base[1]}"
                    if k in lines.keys():
                        lines[k].append(EXTRACTED_FQS_DIR+basename)
                    else:
                        lines[k] =  [EXTRACTED_FQS_DIR+basename]
            return lines

        clean_fq_list = "${extracted_fqs_list}".split("],")

        fq_lines = create_lines(clean_fq_list)
        with open("sorted_manifest.csv", 'w', newline='') as outcsv:
            writer = csv.writer(outcsv)
            writer.writerow(["idx", "sample_id", "reads_1", "reads_2", "tax_id"])
            for idx, sample_tax_id in enumerate(fq_lines.keys()):
                key_pair = sample_tax_id.strip(" ").split(",")
                writer.writerow([idx, key_pair[0], fq_lines[sample_tax_id][0], fq_lines[sample_tax_id][1], key_pair[1]])
        """
}