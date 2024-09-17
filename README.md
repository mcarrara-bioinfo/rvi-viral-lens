[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.1-23aa62.svg)](https://www.nextflow.io/) [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/) ![Nf-test](https://img.shields.io/badge/NFtest-%E2%89%A50.8.4-23aa62.svg?labelColor=0000)

[ADD VIRAL PIPELINE LOGO]

The **VIRAL_PIPELINE** is a Nextflow pipeline developed under the context of the [RVI project]() by [GSU]() and its main goal is to identify the presence of Flu, SARS-CoV-2 and RSV and obtain, if possible, high quality consensus sequences for those virus. For more details, check the [[ADD REFERENCE PAPER]]()

---
**Table of Contents**
- [Pipeline Summary](#pipeline-summary)
- [How to Cite](#how-to-cite)
- [Quick Start](#quick-start)
- [Installation](#installation)
  - [dependencies](#dependencies)
  - [build containers](#build-containers)
- [Usage](#usage)
- [Inputs](#inputs)
- [Parameters](#parameters)
- [Unit Tests](#unit-tests)
- [Pipeline components documentation](#pipeline-components-documentation)
- [Licence](#licence)

---

## Pipeline Summary

The pipeline takes a manifest containing  **fastq pairs file** paths and a **kraken detabase** as inputs (check [Inputs section](#inputs) for more details) and outputs a **classification report**, **consensus sequences** and a **collection of intermediate files** (check [[ADD OUTPUT DOCUMENTATION LINK]]() for more details). Here is an broad overview of the pipeline logic

1. **Sort Reads**: The initial step is sort reads using `kraken2` for each fastq pairs according to the database provided. The classified reads is used as input to `kraken2ref` which will generate one pair of fastq files per taxid found.
   - An option to split big files is provided (check [Parameter section](#parameters)).
   - For a more in depth technical description check [`SORT_READS_BY_REF` workflow documentation](./docs/workflow/SORT_READS_BY_REF.md)

2. **Generate Consensus**: After all samples been classified, all references observed for that samples batch are fetch from the `kraken database` (or an arbitrary fasta file provided by the user). The classfied reads are aligned to their respective references (via `bwa`). The alignment is used as input for `ivar` to obtain a consensus sequence.
   - For a more in depth technical description check [`GENERATE_CONSENSUS` workflow documentation](./docs/workflow/GENERATE_CONSENSUS.md)
3. **Compute QC**: QC metrics are computed via `samtools` and a custom script (`qc.py`)
   - For a more in depth technical description check [`COMPUTE_QC_METRICS` workflow documentation](./docs/workflow/COMPUTE_QC_METRICS.md)

4. **SARS-CoV-2 Subtyping**: SARS-CoV-2 subtyping can be done if present on the sample
   - For a more in depth technical description check [`SCOV2_SUNTYPING` workflow documentation](./docs/workflow/SCOV2_SUBTYPING.md)

![](./docs/assets/metro_map.png)

---


## How to Cite

[ADD CITATION]

---
## Quick Start

Assuming [dependencies](#dependencies) are installed on the system:

1. Setup the pipeline:

```bash
# clone the repo
git clone https://gitlab.internal.sanger.ac.uk/malariagen1/viral/viral_pipeline.git
cd viral_pipeline/
# build containers
cd containers/
sudo singularity build base_container.sif baseContainer.sing
sudo singularity build ivar.sif ivarContainer.sing
sudo singularity build pangolin.sif pangolinContainer.sing
sudo singularity build kraken.sif krakenContainer.sing
sudo singularity build kraken2ref.sif kraken2ref.sing
```

2. Run pipeline

You will need a manifest and a kraken database (check [Inputs section](#inputs) and [Usage section](#usage) for more details)

```bash
PIPELINE_CODES=<path to viral pipeline>
MANIFEST=<path to my manifest>
kraken_db_path=<path to my kraken DB>
PIPELINE_CONTAINERS=<path to my containers dir>

nextflow run ${PIPELINE_CODES}/main.nf --manifest ${MANIFEST} \
    --db_path ${kraken_db_path} \
    --outdir ./output/ \
    --containers_dir ${PIPELINE_CONTAINERS} \
    -with-trace -with-report -with-timeline \
    -profile sanger_local \
    -resume
```

---
## Installation

### dependencies

- [Nextflow](https://www.nextflow.io) (tested on v23.10.0)
- [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) (required to use Singularity Containers)

> We strongly recommend to run the pipeline using the containers recipe provided at `containers/` subdir.

If not using containers, all the software needs to be available at run time. A list of the requirement softwares needed can be found [Software Requirements documentation](./doc/requirements.md)

> TODO: add list of software versions in use

### build containers

Recipes for the Singularity container used on this pipeline are available on this repository at `containers/` dir. To build the containers, run the commands bellow.

```{bash}
cd containers/
sudo singularity build base_container.sif baseContainer.sing
sudo singularity build ivar.sif ivarContainer.sing
sudo singularity build pangolin.sif pangolinContainer.sing
sudo singularity build kraken.sif krakenContainer.sing
sudo singularity build kraken2ref.sif kraken2ref.sing
```

> NOTE: Currently we only support running on Singularity local containers, but we should add Docker container registries soon

## Usage

1. **Generate manifest**
For convenience, a script to generate the manifest is provided on this repo

```{bash}
python write_manifest.py ./path/to/my/fastqs_dir/ -fq1_ext my_r1_ext -fq2_ext my_r2_ext
```

2. **Run pipeline**

```{bash}
nextflow run /path/to/rvi_consensus_gen/main.nf --manifest /path/to/my/manifest.csv \
        --db_path /path/to/my/kraken_db \
        --results_dir outputs/ \
        --containers_dir /path/to/my/containers_dir/ \
        -profile sanger_stantard -resume -with-trace -with-report -with-timeline
```

Optionally is possible to start the pipeline from **GENERATE_CONSENSUS**

```{bash}
nextflow run ${CHECKOUT}/main.nf --entry_point consensus_gen \
        --consensus_mnf sorted_manifest.csv
        --results_dir $LAUNCHDIR/outputs/ \
        --containers_dir /path/to/containers_dir/ \
        -profile sanger_standard-resume -with-trace -with-report -with-timeline
```

## Inputs

- `manifest` : Manifest of input fastq file pairs.
  - Must have `sample_id`,`reads_1` and `reads_2` collumns
- `db_path` : Path of kraken2 database

> **NOTE**: if using the `consensus_gen` entry point, the manifest must containing pair ended reads and genome reference files for each pair. must have `sample_id`, `taxid`, `ref_files`, `reads_1`, `reads_2`

## Parameters

- **General**

  - `containers_dir` [DEFAULT =  `containers/` dir of this repository] : By default, the pipeline relies on Singularity containers and __assumes__ all containers are present on this directory and were named on a specific manner
  - `results_dir` [DEFAULT = `$launchDir/output/`] : set where output files should be published. By default, it will write files to an `output/` dir (if not existent, it will be created) at pipeline launch directory.

- **Pipeline flow controls**

  - `entry_point` [Default = `"sort_reads"`] : Defines the entry point of the pipeline. The `sort_reads`
  - `scv2_keyword` [Default = `"Severe acute respiratory syndrome coronavirus 2"`] : keyword to identified samples with SARS-CoV-2. This string is obtained from the kraken2 database, therefore, it should be in line with the database in use.
  - `do_scov2_subtyping` [Default = `true`] : Switch SARS-CoV-2 on and off.

- **Specific processes**

  - `get_taxid_references`:
    - `db_library_fa_path` [Default = `null`]: A fasta file containing reference sequences for the taxids present on the kraken database. It assumes there is a `$db_`
  - `raken2ref`:
    - `k2r_fq_load_mode` [Default = `"full"`] : kraken2ref load fastq into memory mode ["full"/"chunk"]
    - `k2r_max_total_reads_per_fq` [Default = `10000000`] : set maximum number of reads to be accepted any classified fastq file pair. Files excedding that limit will be splitted before `run_k2r_dump_fastq`.
    - `k2r_dump_fq_mem` [Default = `"6 GB"`] : Memory to be requested by `run_k2r_dump_fq`, adjust according to `k2r_max_total_reads_per_fq`.
    - `min_reads_for_taxid` [Default = `100`] : min reads number to be considered.
  - `ivar`
    - `ivar_min_depth` [Default = `10`] : Minimum depth to call consensus
    - `ivar_freq_threshold` [Default = `0.75`] : Minimum frequency threshold(0 - 1) to call consensus.

> NOTE: `Kraken2ref` have a escalation memory strategy based on linear regression, check [k2r_memory_escalation documentation](./docs/k2r_memory_escalation.md) for more details.

> TODO: check if all pipeline parameters are covered and add missing ones

## Unit Tests
---
The workflow & process unit tests for this pipeline are written in the [nf-test](https://www.nf-test.com/) (`v0.8.4`) Nextflow testing framework.

**Running Tests**

The following command if entered from the repository top-level directory can be used to execute all of the per-process & per-workflow unit tests:

**Run all tests**

```{bash}
nf-test test ./
```

**Run individual module/workflow test**

```{bash}
nf-test test tests/<modules or workflows>/<module_to_test>.nf.test
```

---
## Pipeline components documentation

For a more in depth technical documentation of all the processes and workflows can be found at `docs/` dir, here you will find documentation for:

**Processes**

- [bwa_alignment.nf](./docs/modules/bwa_alignment.md)
- [get_taxid_references.nf](./docs/modules/get_taxid_references.md)
- [run_ivar.nf](./docs/modules/run_ivar.md)
- [run_kraken.nf](./docs/modules/run_kraken.md)
- [run_kraken2ref_and_pre_report.nf](./docs/modules/run_kraken2ref_and_pre_report.md)
- [run_pangolin.nf](./docs/modules/run_pangolin.md)
- [run_qc_script](./docs/modules/run_qc_script.md)
- [write_classification_report](./docs/modules/write_classification_report.md)

**Workflow**

- [SORT_READS_BY_REF.nf](./docs/workflow/SORT_READS_BY_REF.md)
- [GENERATE_CONSENSUS.nf](./docs/workflow/GENERATE_CONSENSUS.md)
- [COMPUTE_QC_METRICS.nf](./docs/workflow/COMPUTE_QC_METRICS.md)
- [SCOV2_SUBTYPING.nf](./docs/workflow/SCOV2_SUBTYPING.md)

---

## Licence

[ADD LICENSE]
