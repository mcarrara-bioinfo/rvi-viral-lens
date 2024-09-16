# Configuration Documentation: `nextflow.config`

## Overview

This is the main Nextflow configuration file is used to define workflow parameters, resource handling, containerization options, and process execution profiles. It includes configurations for running specific tasks such as sorting reads, generating consensus sequences, and subtyping SARS-CoV-2. The configuration sets Singularity containerization and compute resources across different execution environments (`local` and `LSF`).

It includes two other config files to handle [singularity containers] from `./conf/containers.nf` and resource handling for kraken2ref `./conf/resource_handling.nf`

## Key Components

1. **Parameters (`params`)**: Defines workflow entry points, database paths, and task-specific settings.
2. **Containerization**): Determines the use of Singularity for process execution.
3. **Process Settings (`process`)**: Defines caching and execution behavior.
4. **Profiles**: Configures specific environments and resource requirements, including containerization and job execution details for running on different platforms. Currently, only a profile for execution under Wellcome Sanger Instite computational environment.

---

## Configuration Sections

### 1. **Parameters (params)**

The `params` block defines key user-modifiable settings for the workflow.

- `entry_point`: Entry point for the workflow. It can be set to `"sort_reads"` or `"consensus_gen"`.
  - Default: `"sort_reads"`.
- `manifest`: Path to the manifest file (default: `null`).

**Kraken Database Parameters**

- `db_path`: Path to the Kraken database.
- `db_library_fa_path` (OPTIONAL): Path to the Kraken database library FASTA file.
  - By default, it assumes there is a `${params.db_path}/library/library.fna`.

**Kraken2Ref Handling**

- `k2r_fq_load_mode`: Loading mode for Kraken2 fastq files (either `full` or `chunks`).
  - Default: `"full"`.

- `k2r_max_total_reads_per_fq`: Maximum number of reads to process per fastq file.
  - Default: `10,000,000`.
- `k2r_dump_fq_mem`: Memory allocated for dumping fastq files.
  - Default: `"6 GB"`.

**Kraken2ref Report Filter**

- `min_reads_for_taxid`: Minimum number of reads required to assign a taxonomic ID.
  - Default: `100`.

**iVar Parameters**

- `ivar_min_depth`: Minimum depth for consensus calling.
  - Default: `10`.
- `ivar_freq_threshold`: Sets base frequency threshold for consensus calling.
  - Default: `0.75`.

**Virus Subtyping**

- `scv2_keyword`: Keyword to identify SARS-CoV-2 sequences. Any taxid name equal to the string set by this parameter will be considered as SCOV2 and subjected to specific SARS-CoV-2 subtyping.
  - Default: `"Severe acute respiratory syndrome coronavirus 2"`.

- `do_scov2_subtyping`: Boolean flag to enable or disable SARS-CoV-2 subtyping via Pangolin. Check [SARS0COV2 documentation](./workflow/SCOV2_SUBTYPING.md) for more details
  - Default: `true`.

**Consensus Generation**

- `consensus_mnf` [OPTIONAL]: Path to the consensus manifest file (default: `null`).

**Output Directory**

- `results_dir`: Directory where output files will be published.
  - Default: `"$launchDir/output/"`.

### 2. **Containerization**

Currently, the pipeline only provide Singularity containers. Check the [specific containers documentation](./containers.md).

**Docker**

- `enabled`: Flag to enable or disable Docker.
  - Default: `false`.

**Singularity**

- `enabled`: Flag to enable or disable Singularity.
  - Default: `true`.

By default, the current version of the pipeline assumes all singularity containers are available on specific paths with specific names defined at `./conf/containers.config`. 
- Currently, `"$projectDir/containers/"` is the default location for the containers, it can be changed by the user using `containers_dir` parameter.


#### Containers
All containers used by this pipeline recipes can be found at `./containers/` dir of this repository. The current containers in use are:

- `base_container.sif`: the default container for all processes unless overridden.
- `ivar.sif`: container to be used for processes using the `ivar` label
- `kraken.sif`: container to be used for processes using the `kraken` label
- `pangoling.sif`: container to be used for processes using the  `pangolin` label
- `kraken2ref.sif`: container to be used for processes using the  `kraken2ref` label.

### 3. Process Settings (`process`)

- `cache='lenient'`: Defines the cache behavior for processes, allowing cached results to be reused even when minor changes occur.
- `executor='local'`: The default executor is set to local, meaning processes will run on the local machine unless otherwise specified.

### 4. Profiles
Profiles define environment-specific configurations. Currently, there is a single predefined profile: `sanger_standard`.

#### `sanger_standard` Profile

This profile is optimized for running jobs on the Sanger Institute's infrastructure and includes settings for using Singularity, job execution on the LSF scheduler, and handling of specific processes.

**Singularity Settings**

- `autoMounts`: Automatically mounts paths required for job execution.
- `cacheDir`: Cache directory set to "`$PWD`" (the current working directory).
- `runOptions`: Singularity run options to bind necessary paths (`/lustre`, `/nfs`, `/software`, `/data/`).

**Process-Specific Settings**
Inherited from the global process configuration
- `cache='lenient'`: this makes `resume` more tolerant to changes in files attributes, such as timestamps. This is handy when using distributed file system which uses [Network File System protocols](https://en.wikipedia.org/wiki/Network_File_System), check [Nextflow documentation](https://www.nextflow.io/docs/latest/cache-and-resume.html#inconsistent-file-attributes) for more details.
- `executor='local'`: Default executor is local unless overridden.

**LSF Job Execution**
Certain processes are configured to run on the LSF scheduler with specific resource allocations:

- `run_k2r_sort_reads`:
  - Executor: `lsf`
  - CPU cores: `1`

- `run_k2r_dump_fastqs_and_pre_report`:
  - Executor: `lsf`
  - CPU cores: `1`

- `run_kraken`:
  - Executor: `lsf`
  - CPU cores: `16`
  - Memory: `1 GB`

**Job Naming and Memory Management**

- `jobName`: Custom job name format for LSF jobs, based on the task name and tag (`"RVI-preprocess - $task.name - $task.tag"`).

- `perJobMemLimit=true`: Ensures that memory limits are set per job.
