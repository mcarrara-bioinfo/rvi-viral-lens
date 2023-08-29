#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
include {load_json_file} from './workflow/load_json.nf'

// Main entry-point workflow
workflow {

    // === 1 - Process input ===
    // 1.1 - Get reads files 
    reads_channel_fqgz = channel.fromFilePairs("${params.reads_dir}/*_R{1,2}*.fq.gz")
    reads_channel_fagz = channel.fromFilePairs("${params.reads_dir}/*_R{1,2}*.fastq.gz")
    reads_ch = reads_channel_fagz.concat(reads_channel_fqgz)
    reads_ch.view()
    // 1.2 - load virus settings
    virus_resources = load_json_file(params.virus_resources_json)
    //println(virus_resources.sarscov2.reference_genome[0])

    // ==========================

    // Do whole reference genome mapping
    //bwa-mem
    //samtools-sort
    //samtools mpileup
    
    // Do genome consensus generation
    //ivar - > consensus
    
    // Do alignment of consensus to genome to ref

    // Do species specific

    // general

    // SARS-CoV-2

    // Flu


}
