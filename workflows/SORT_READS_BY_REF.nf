include {run_kraken} from '../modules/run_kraken.nf'
include {get_taxid_reference_files} from '../modules/get_taxid_references.nf'
include {run_kraken2ref_and_pre_report} from '../modules/run_kraken2ref_and_pre_report.nf'

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

        // 0 - create channel from input per-sample manifest
        mnf_ch = parse_clean_mnf_meta(clean_mnf)
        // -------------------------------------------------//


        // 1 - run kraken and get outputs
        run_kraken(mnf_ch, params.db_path)
        // -------------------------------------------------//
        
        // 2 - obtain unique taxid reference fasta file and metadata from kraken db 

        // check if library fasta was set, 
        // if not, assume is at library/library.fna on the kraken db dir 
        if (params.db_library_fa_path == null){
            library_fa_path = "${params.db_path}/library/library.fna"
            log.warn("No db_library_fa_path set, assuming ${library_fa_path} exists")
        } else {
            library_fa_path = params.db_library_fa_path
        }
        // -------------------------------------------------//

        // 3 - reads per taxid on fastq file

        // drop unclassified fq filepair
        run_kraken.out // tuple (meta, kraken_output, [classified_fq_filepair], [unclassified_fq_filepair], kraken_report)
            | map {it -> tuple(it[0], it[1], it[2], it[4])} // tuple (meta, kraken_output, [classified_fq_filepair], kraken_report)
            | set {sort_reads_in_ch}

        // -------------------------------------------------//

        // 4 - run kraken2ref and collect all per-sample per-taxon fq filepairs
        run_kraken2ref_and_pre_report(sort_reads_in_ch)
        sample_pre_report_ch = run_kraken2ref_and_pre_report.out.report_file
        
        // prepare channel to be emitted
        run_kraken2ref_and_pre_report.out.fq_files // tuple (meta, [id_taxid_R{1,2}.fq])
            | map {meta, reads ->
                // group pairs of fastqs based on file names, and add new info to meta
                reads
                    .groupBy { filePath -> filePath.getName().tokenize("_")[0..-2].join(".")}
                    .collect { identifier, paths ->[[identifier, paths]]}
            // this map returns a channel with a single value,
            // which is a list with all the ids and files.
            }
            | flatten() // flat the list [id_a, fqa1, fqa2, idb, fqb1, fqb2, ...]
            | collate(3) // tuple (id.taxid, fq_1, fq_2)
            | map { it ->
                // rebuild meta and reads structure
                meta = [sample_id:it[0].tokenize(".")[0..-2].join("_"), //run.lane.sample_info
                        taxid:it[0].tokenize(".")[-1]] //taxid
                meta.id = "${meta.sample_id}.${meta.taxid}"
                reads = [it[1], it[2]]
                [meta, reads]
            }
            | set {pre_sample_taxid_ch}

        // 5 - get references files
        pre_sample_taxid_ch
            .map{meta, reads -> tuple(meta.taxid)}
            .collect()
            .flatten()
            .unique()
            .set {unq_taxid_ch}

        get_taxid_reference_files(unq_taxid_ch, library_fa_path)
        get_taxid_reference_files.out// tuple (taxid, [ref_files])
            .set{taxid_ref_files_map_ch}
        
        // add reference files to meta
        pre_sample_taxid_ch
            | map {meta, reads ->
                [meta.taxid, meta, reads]
            }
            | combine(taxid_ref_files_map_ch, by:0) // [taxid, meta, reads, ref_files]
            | map {taxid, meta, reads, ref_files ->
                meta.ref_files = ref_files
                [meta, reads]
            }
            | set {sample_taxid_ch}

    emit:
        sample_taxid_ch // tuple (meta, reads)
        sample_pre_report_ch // tuple (meta,)
}

def check_sort_reads_params(){
    def errors = 0
    // was the kraken database provided?
    if (params.db_path == null){
        log.error("No kraken database path provided")
        errors +=1
    }

    // TODO: if provided, check if it is a valid path

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