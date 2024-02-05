process compute_number_of_aligned_reads {
    tag "${meta.id}"
    input:
        tuple val(meta), path(aligned_bam)
    output:
        tuple val(meta), stdout

    script:
$/
#!/usr/bin/python3
from io import StringIO
import pysam

def read_tsv_as_dict_from_stringio(stringio_object):
    data_dict = {'chromosome': [], 'length': [], 'mapped': [], 'unmapped': []}

    for line in stringio_object.getvalue().splitlines():
        # Split each line into fields using tab as the delimiter
        fields = line.strip().split('\t')[0:4]
        # skip line with chr as "*"
        if fields[0] == "*":
            continue

        # Update the dictionary with values from each column
        for key, value in zip(data_dict.keys(), fields):
            data_dict[key] = value #.append(value.replace(',','.'))

    return data_dict

#run pysam idxstats
idx = StringIO(pysam.idxstats("${aligned_bam}"))
tsv_data_dict = read_tsv_as_dict_from_stringio(idx)
print(f"{tsv_data_dict['mapped']}|{tsv_data_dict['unmapped']}|")
/$
}