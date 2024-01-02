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

def parse_clean_mnf_meta(consensus_mnf) {

    def mnf_ch = Channel.fromPath(consensus_mnf)
                    | splitCsv(header: true, sep: ',')
                    | map {row -> 
                        // set meta
                        meta = [id: row.sample_id]
                        // set files
                        reads = [row.reads_1, row.reads_2]
                        // declare channel shape
                        [meta, reads]
                    }
    return mnf_ch // tuple(sample_id, [fastq_pairs])
}

workflow SORT_READS_BY_REF {
    take:

        clean_mnf // path to clean manifest

    main:
        // create channel from input per-sample manifest
        mnf_ch = parse_clean_mnf_meta(clean_mnf)

        // run kraken and get outputs
        run_kraken(mnf_ch, params.db_path, params.results_dir)

        // drop unclassified fq filepair
        run_kraken.out // tuple (meta, kraken_output, [classified_fq_filepair], [unclassified_fq_filepair], kraken_report)
            | map {it -> tuple(it[0], it[1], it[2], it[4])} // tuple (meta, kraken_output, [classified_fq_filepair], kraken_report)
            | set {sort_reads_in_ch}

        // run krakentools and collect all per-sample per-taxon fq filepairs
        sort_reads_with_krakentools(sort_reads_in_ch)

        //sort_reads_with_krakentools.out.set {sample_taxid_ch}


        sort_reads_with_krakentools.out // tuple (meta, [id.tax_id.extracted{1,2}.fq])
            | map {meta, reads ->
                // group pairs of fastqs based on file names, and add new info to meta
                reads
                    .groupBy { filePath -> filePath.getName().tokenize(".")[0..1].join(".")}
                    .collect { identifier, paths ->[[identifier, paths]]}
            // this map returns a channel with a single value,
            // which is a list with all the files.
            }
            | flatten() // flat the list [id_a, fqa1, fqa2, idb, fqb1, fqb2, ...]
            | collate(3) // tuple (id.taxid, fq_1, fq_2)
            | map { it ->
                // rebuild meta and reads structure
                meta = [sample_id:it[0].tokenize(".")[0],
                        taxid:it[0].tokenize(".")[1]]
                meta.id = "${meta.sample_id}.${meta.taxid}"
                reads = [it[1], it[2]]
                [meta, reads]
            }
            | set {sample_taxid_ch}

        // ---------------------------- //


        // this should still be used, but should not be go downstream
        //sort_reads_with_krakentools.out.collect().view()
        //output_mnf_ch = sort_reads_with_krakentools.out.map{it -> it[1]}.collect()
        
        //output_mnf_ch.view()
        // write a per-sample per-taxon manifest
        //write_sorted_manifest(output_mnf_ch)
        //consensus_mnf = write_sorted_manifest.out
    
    emit:
        sample_taxid_ch // tuple (meta, reads) 
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