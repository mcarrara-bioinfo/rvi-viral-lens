# Changelog

All notable changes to this project will be documented in this file.
## Unreleased

### Changed

- **[improvement]**: Mpileup output retained by run_ivar & used by the QC script for calculating % genome coverage.
- **[improvement]**: Removed unnecessary code from qc.py and run_qc.nf including the plot generation.
- **[improvement]**: Modified qc.py to read input files from command line including using samtools flagstat for read counts.
- **[improvement]**: Unit test files and snapshot files for run_ivar, run_qc_script, and GENERATE_CONSENSUS to account for changes

### Added

- **[improvement]**: Update container of kraken2ref from v2.0.0 to v2.05.0
- **[added]**: add a paramter to set the pooling method for kraken to ref (default kmeans)
- **[added]**: Container for the run_qc process
- **[added]**: Unit test for COMPUTE_QC_METRICS
- **[added]**: Mpileup test data
- **[improvement]**: implement k2r release new features
- **[improvement]**: split fastq files if higher than a set numbers of reads per fq

## [0.2.2] - 2024-08-02

### Changed

- **[improvement]**: Ivar module now follows Artic pipeline standards

### Added

- **[added]**: add LSF memory escalation strategy for kraken2ref 
- **[added]**: Columns Virus_Taxon_ID, Virus, Species, Reference_Taxon_ID, Selected_Reference added/populated to classification report

## [0.2.1] - 2024-06-20

### Fixed

- **[bug]**: Classification report generation would crash if ' was present in output report file lines
- **[bug]**: Independent workflow stanza for GENERATE_CLASSIFICATION_REPORT.nf was outdated / broken


## [0.2.0] - 2024-05-29

### Fixed

- **[bug]**: Classification report and pre report parsing errors fixed

### Changed

- **[improvement]**: Remove mpileup repeated command calls on ivar process.
- **[improvement]**: Remove redundant processes, rewiring and tiding up code base.
- **[improvement]**: Qc metrics using the same method of the Artic pipeline
- **[improvement]**: add kraken2ref as the new reference selection tool
- **[updated]**: Unit tests adapted to new channel and processes structure

### Added

- **[added]**: A script (`k2r_report.py`) was added to generate a pre report file from k2r software
- **[added]**: unit tests for new `run_kraken2ref_and_pre_report.nf`

## [0.1.0] - 2024-02-05

### Added

- **[added]**: Viral subtyping and classification reports routines integrated to pipeline
- **[added]**: `Percentage Coverage` and `number of mapped reads` are now computed at a new QC metrics workflow
- **[added]**: new workflow (`SUBTYPE_AND_SEGMENT_FLU.nf`) attempts to retrieve the flu subtype and segment from kraken report file and populates the meta with these values
- **[added]**: new module (`retrieve_flu_subtype_and_segment.nf`) attempts to parse out the flu subtype and segment from the kraken report file and sets these values to Null if nothing retrieved
- **[added]**: QC metrics workflow, currently computes reads depth and percentage genome coverage
- **[added]**: SARS-CoV-2 sequences subtyping via pangolin
- **[added]**: branching `GENERATE_CONSENSUS` workflow output for viral subtyping routines
- **[added]**: new parameter (`min_reads_for_taxid`)to set a treshold for minimum number of reads assigned for a taxid to be considered
- **[added]**: new workflow and module (`GENERATE_CLASSIFICATION_REPORT and write_classification_report`) to generate a classifcation report file
- **[added]**: unit tests written in `nf-test` covering modules, workflows and pipeline

### Changed

- **[improvement]**: `taxid` respective `rank` and `name` are available on meta
- **[improvement]**: Taxid reference fasta files for consensus sequence are obtained from kraken database
- **[improvement]**: Channels now rely on Meta Mapping
- **[improvement]**: Output folder now have the following structure `output/<sample_id>/<taxid>`
- **[improvement]**: `write_manifest.py` relies on glob expression

### Fixed

- **[bug]**: Samples with no taxids to process raises a warning and are now ignored

### Removed

- **[Removed]**: writing manifest process removed from `SORT_READS_BY_REF.nf` 
- **[Removed]**: json resource files and fasta files provided on the repo

---
## [0.0.1] - 2023-12-01

This is the first prototype versioning. This pipeline 1) sort reads via Kraken and 2) generate consensus sequences using ivar.
