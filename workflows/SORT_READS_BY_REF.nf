include {run_kraken} from '../modules/run_kraken.nf'
include {sort_reads_with_krakentools} from '../modules/sort_reads.nf'
include {get_taxid_reference_files} from '../modules/get_taxid_references.nf'
include {filter_taxids} from '../modules/filter_taxid.nf'

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


        // 2 - collect all taxid observed
        run_kraken.out
            | map {it -> tuple(it[0], it[-1])}
            | set {kraken_reports_ch}

        filter_taxids(kraken_reports_ch)

        filter_taxids.out // [meta, taxid_string]
            .map {meta, taxid_str -> taxid_str.tokenize(" ")}
            .collect()
            .flatten()
            .unique()
            .set {unq_taxid_ch}
        // -------------------------------------------------//

        // 3 - obtain unique taxid reference fasta file from kraken db 

        // check if library fasta was set, 
        // if not, assume is at library/library.fna on the kraken db dir 
        if (params.db_library_fa_path == null){
            library_fa_path = "${params.db_path}/library/library.fna"
            log.warn("No db_library_fa_path set, assuming ${library_fa_path} exists")
        } else {
            library_fa_path = params.db_library_fa_path
        }

        get_taxid_reference_files(unq_taxid_ch, library_fa_path)

        get_taxid_reference_files.out// tuple (taxid, [ref_files])
            .set{taxid_ref_files_map_ch}
        // -------------------------------------------------//

        // 4 - reads per taxid on fastq file 
        
        // drop unclassified fq filepair
        run_kraken.out // tuple (meta, kraken_output, [classified_fq_filepair], [unclassified_fq_filepair], kraken_report)
            | map {it -> tuple(it[0], it[1], it[2], it[4])} // tuple (meta, kraken_output, [classified_fq_filepair], kraken_report)
            | set {sort_reads_in_ch}
        
        // run krakentools and collect all per-sample per-taxon fq filepairs
        sort_reads_with_krakentools(sort_reads_in_ch)

        sort_reads_with_krakentools.out.set {sample_taxid_ch}
        // -------------------------------------------------//


        // 5 - prepare channel to be emitted
        sort_reads_with_krakentools.out // tuple (meta, [id.tax_id.extracted{1,2}.fq])
            | map {meta, reads ->
                // group pairs of fastqs based on file names, and add new info to meta
                reads
                    .groupBy { filePath -> filePath.getName().tokenize(".")[0..1].join(".")}
                    .collect { identifier, paths ->[[identifier, paths]]}
            // this map returns a channel with a single value,
            // which is a list with all the ids and files.
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
            // add reference files to meta
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