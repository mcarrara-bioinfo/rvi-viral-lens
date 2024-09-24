process get_taxid_reference_files{
    /*
    * Fetch Fasta Sequence Files for a Given Taxid
    *
    * The get_taxid_reference_files process is designed to extract reference 
    * sequences corresponding to a specified taxonomic ID (taxid) from a 
    * larger FASTA file. This process retrieves the relevant sequence, 
    * writes it to an output file, and indexes it using BWA (Burrows-Wheeler 
    * Aligner). This step is essential for downstream analysis that requires
    * taxon-specific reference sequences. In the context of the pipeline, it
    * is used to extract all the reference files from the Kraken database 
    * which were observed on that input batch.
    *
    * check docs/modules/get_taxid_references.md for more extensive documentation
    */

    tag "${taxid}"
    publishDir "${params.results_dir}/reference_files/${taxid}/", mode: 'copy'


    input:
        val(taxid)
        val(kraken_db_library_path)
    output:
        tuple val(taxid), path("${taxid}.fa*"), optional : true

    script:
$/
#!/usr/bin/python3

import subprocess 

taxid_to_find = "${taxid}"
output_file = "${taxid}.fa"
source_fna_path = "${kraken_db_library_path}"

with open(source_fna_path, "r") as source_file:
    header = ""
    seq = ""
    found = False
    nfound = 0
    for line in source_file:
        line = line.strip()
        if line.startswith(">"):
            if found:
                break  # Exit loop after finding the first sequence matching taxid
            header = line
            if f"|{taxid_to_find}|" in header:
                nfound +=1
                found = True
        elif found:
            seq += line

# Write output only if a matching sequence is found
print(f"{nfound} sequences for {taxid_to_find}")
if found:
    with open(output_file, "w") as output:
        output.write(header + "\n" + seq + "\n")

    # Run bwa index on the output file
    subprocess.run(["bwa", "index", output_file])
/$
}