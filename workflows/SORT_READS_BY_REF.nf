include {run_kraken} from '../modules/run_kraken.nf'
include {sort_reads_with_krakentools} from '../modules/sort_reads.nf'
include {write_sorted_manifest} from '../modules/write_sorted_manifest.nf'

//Antonio's code (https://gitlab.internal.sanger.ac.uk/malariagen1/rvi/rvi_consensus_gen/-/blob/main/workflow/GENERATE_CONSENSUS.nf?ref_type=heads) (slightly edited)
def parse_clean_mnf(consensus_mnf) {

    def mnf_ch = Channel.fromPath(consensus_mnf, checkIfExists: true)
                        | splitCsv(header: true, sep: ',')
                        | map { row ->
                            tuple(
                                row.sample_id,
                                [row.reads_1, row.reads_2]
                            )
                        }
    return mnf_ch // tuple(run_id, sample_id, [fastq_pairs])
}

workflow SORT_READS_BY_REF {
    take:

        clean_mnf // path to clean manifest

    main:
        // create channel from input per-sample manifest
        mnf_ch = parse_clean_mnf(clean_mnf)

        // run kraken and get outputs
        run_kraken(mnf_ch, params.db_path, params.results_dir)
        kraken_out_ch = run_kraken.out // tuple (sample_id, kraken_output, [classified_fq_filepair], [unclassified_fq_filepair], kraken_report)

        // drop unclassified fq filepair
        kraken_out_ch
        | map {it -> tuple(it[0], it[1], it[2], it[4])} // tuple (sample_id, kraken_output, [classified_fq_filepair], kraken_report)
        | set {sort_reads_in_ch}

        // run krakentools and collect all per-sample per-taxon fq filepairs
        sort_reads_with_krakentools(sort_reads_in_ch)
        output_mnf_ch = sort_reads_with_krakentools.out.collect()

        // write a per-sample per-taxon manifest
        write_sorted_manifest(output_mnf_ch)
        consensus_mnf = write_sorted_manifest.out

    emit:
        consensus_mnf // path to /work dir copy of output manifest

}

def check_sort_reads_params(){
    def errors = 0
    // was the manifest provided?
    if (params.manifest == null){
        log.error("No manifest provided")
        errors +=1
    }
    // if yes, is it a file which exists? 
    if (params.manifest){
        manifest_file = file(params.manifest)
        if (!manifest_file.exists()){
            log.error("The manifest provided (${params.manifest}) does not exist.")
            errors += 1
        }
        //TODO
        //else {
        //    validate_manifest(params.manifest)
        //}
    }
    return errors
}