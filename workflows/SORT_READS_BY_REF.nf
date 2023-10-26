include {run_kraken} from '../modules/run_kraken.nf'
include {sort_reads_with_krakentools} from '../modules/sort_reads.nf'
include {write_sorted_manifest} from '../modules/write_sorted_manifest.nf'

//Antonio's code (https://gitlab.internal.sanger.ac.uk/malariagen1/rvi/rvi_consensus_gen/-/blob/main/workflow/GENERATE_CONSENSUS.nf?ref_type=heads) (slightly edited)
def parse_clean_mnf(consensus_mnf) {

    def mnf_ch = Channel.fromPath(consensus_mnf, checkIfExists: true)
                        | splitCsv(header: true, sep: ',')
                        | map { row ->
                            tuple(
                                row.run_id,
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

        mnf_ch = parse_clean_mnf(params.manifest)
        mnf_ch.view()


        run_kraken(mnf_ch, params.db_path, params.outdir)
        kraken_out_ch = run_kraken.out

        // kraken_out_ch.view()

        kraken_out_ch
        | map {it -> tuple(it[0], it[1], it[2], it[4])} //tuple (sample_id, kraken_output, [classified_fq_filepair], kraken_report)
        | set {sort_reads_in_ch}

        // sort_reads_in_ch.view()

        sort_reads_with_krakentools(sort_reads_in_ch)

        output_mnf_ch = sort_reads_with_krakentools.out
        // output_mnf_ch.view()

        write_sorted_manifest(output_mnf_ch)
        consensus_mnf = write_sorted_manifest.out

    emit:
        consensus_mnf

}
