#!/usr/bin/env python3
import re
import json
import pandas as pd
import argparse

def get_args():
    parser = argparse.ArgumentParser("Script to convert kraken2ref JSON to TSV report.")
    parser.add_argument("-i", "--input_json", type = str, required = True, help = "The JSON file produced by kraken2ref.")
    parser.add_argument("-r", "--report", type = str, required = True, help = "The kraken2 taxonomic report.")
    parser.add_argument("-out_suffix", type = str, required=False, default = ".viral_pipe.report.tsv", 
        help="report filname suffix (default: '.viral_pipe.report.tsv'")
    return parser

def regex_subtyping(regex, search_string):
    pattern = re.compile(regex)
    found = pattern.search(search_string)
    if found:
        return found.group()
    else:
        return None


def get_report(in_file, report_file, out_suffix=".viral_pipe.report.tsv"):

    ## read in kraken report and kraken2ref JSON
    report_df = pd.read_csv(report_file, sep = "\t", header = None)
    ref_json = json.load(open(in_file))

    ## collect tax_ids for selected reference taxa
    selected_ref_taxa = ref_json["metadata"]["selected"]

    ## collect tax_ids for Species associated with each chosen reference tax_id
    virus_ids = [int(ref_json["outputs"][k]["source_taxid"]) for k in ref_json["outputs"].keys()]

    ## collect the human-readable names associated with the viru_ids and selected_ref_taxa
    ## from kraken report
    virus_names_df = report_df[report_df[4].isin(virus_ids)]
    virus_desc_names_dict = dict(zip(virus_names_df[4], [i.strip() for i in virus_names_df[5]]))

    subset_df = report_df[report_df[4].isin(selected_ref_taxa)]
    desc_names_dict = dict(zip(subset_df[4], [i.strip() for i in subset_df[5]]))

    ## get the sample_id
    sample_id = ref_json["metadata"]["sample"]

    ## initialise output data structure
    report_output = {}
    for idx, selected_ref in enumerate(selected_ref_taxa):
        ## get data specific to the selected_ref
        selected_data_dict = ref_json["outputs"][str(selected_ref)]

        ## identify its Species and get Species name
        virus = selected_data_dict["source_taxid"]
        virus_name = None
        virus_name = virus_desc_names_dict[virus]

        ## get the name of the actual reference that was chosen
        ref_name = desc_names_dict[selected_ref]

        ## if selected_ref is ANY flu segment, collect its subtype info
        ## THIS IS THE SUBTYPE OF THE CHOSEN REFERENCE ONLY
        ## THIS MAY NOT BE THE SAME AS THE SUBTYOE OF THE SAMPLE ITSELF
        generic_subtype = regex_subtyping("H[0-9]+N[0-9]+", ref_name)

        ## if selected_ref is Flu Segment 4 or 6, collect subtype info
        ## THIS WILL BE THE SUBTYPE OF THE SAMPLE
        subtype = "None"
        segment = "None"
        if generic_subtype is not None:
            segment = regex_subtyping("(?<=segment )[0-9]", ref_name)
        if segment == 4:
            subtype = regex_subtyping("H[0-9]+", ref_name)
        if segment == 6:
            subtype = regex_subtyping("N[0-9]+", ref_name)

        ## collect number of reads written to each fastq filepair
        if "per_taxon" in ref_json["metadata"]["summary"].keys():
            num_reads = ref_json["metadata"]["summary"]["per_taxon"][str(selected_ref)]
        else:
            num_reads = None

        if "Alphainfluenzavirus" in virus_name:
            report_name = "Influenza A Virus"
        elif "Betainfluenzavirus" in virus_name:
            report_name = "Influenza B Virus"
        else:
            report_name = virus_name

        ## populate output data dict
        report_output[idx] = {
                                "sample_id": sample_id,
                                "virus": virus,
                                "virus_name": virus_name,
                                "selected_taxid": selected_ref,
                                "ref_selected": ref_name,
                                "sample_subtype": subtype,
                                "flu_segment": segment,
                                "virus_subtype": generic_subtype,
                                "parent_selected": selected_data_dict["parent_selected"],
                                "num_reads": num_reads,
                                "report_name": report_name
                            }

    ## once iteration over all chosen refs is done
    ## convert dict of dicts to tabular format and write to tsv
    output_df = pd.DataFrame.from_dict(report_output.values())
    # if not empty, sort dataframe and write csv
    if len(output_df) > 0:
        output_df.sort_values("selected_taxid").to_csv(f"{sample_id}{out_suffix}", sep = "\t", header=True, index = False)
    # if empty, write a csv file
    else:
        output_df.to_csv(f"{sample_id}{out_suffix}", sep = "\t", header=True, index = False)
def main():
    args = get_args().parse_args()
    get_report(args.input_json, args.report, out_suffix=args.out_suffix)

if __name__ == "__main__":
    main()
