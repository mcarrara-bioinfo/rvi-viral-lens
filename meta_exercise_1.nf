#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

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
clean_mnf = "/lustre/scratch126/gsu/team112/personal/ad45/rvi_box/DT1864_meta_pattern_workshop/input_fq/manifest.csv"

workflow {
    // --- META MUST BE INSERTED HERE --- //
    mnf_ch = parse_clean_mnf(clean_mnf)
    mnf_ch.view()
    //-----------------------------------//

}