# Workflow Documentation: `SCOV2_SUBTYPING`

## Overview

The `SCOV2_SUBTYPING` workflow is designed to determine the SARS-CoV-2 lineage (subtype) of consensus sequences using the [PANGOLIN](https://github.com/cov-lineages/pangolin) tool. This workflow takes in a channel of consensus sequences along with their metadata, runs the PANGOLIN lineage classification, and outputs updated metadata with the assigned lineage.

## Key Processes

- **SARS-CoV-2 Lineage Classification**: The workflow uses the PANGOLIN tool to classify each consensus sequence into a SARS-CoV-2 lineage.
- **Metadata Update**: The assigned lineage is added to the sample's metadata.

## Workflow Inputs and Outputs

### Inputs

- **Consensus Sequence Channel (`consensus_seq_ch`)**: A channel containing tuples of metadata and consensus sequences. Each tuple should include:
  - `meta`: Metadata for each sample, which must include the following keys:
    - `id`: Unique identifier for the sample.
    - `taxid`: Taxonomic ID of the sample.
    - `sample_id`: Sample identifier.
  - `consensus_seq`: The consensus sequence (FASTA file) of the sample.

### Outputs

- **SCOV2 Subtype Channel (`scov2_subtype_out_ch`)**: A channel emitting updated metadata for each sample, now including the assigned SARS-CoV-2 lineage (virus subtype).

## Workflow Steps

1. **SARS-CoV-2 Lineage Classification**
The workflow starts by classifying the consensus sequences into SARS-CoV-2 lineages:

- The run_pangolin process is executed on the `consensus_seq_ch` channel.
- Each consensus sequence is processed by the PANGOLIN tool, which assigns a lineage based on the sequence data.

2. **Update Metadata with Lineage**
After the PANGOLIN classification:

- The lineage information (`virus_subtype`) is added to the sample's metadata (`meta`).
- The updated metadata, including the lineage, is output as a `tuple (tuple(meta))` and collected into the `scov2_subtype_out_ch` channel.

3. **Emit Output**
The final output channel (`scov2_subtype_out_ch`) emits the updated metadata with the SARS-CoV-2 lineage classification.

## Example Input Channel
An example input channel (consensus_seq_ch) might look like this:

```groovy
// Example of creating the consensus_seq_ch channel
consensus_seq_ch = Channel.of(
    [ [id: 'sample_001', taxid: '2697049', sample_id: 'sample_001'], 'consensus_sample_001.fasta' ],
    [ [id: 'sample_002', taxid: '2697049', sample_id: 'sample_002', 'consensus_sample_002.fasta' ]]
)
```
