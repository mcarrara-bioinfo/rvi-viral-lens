# Changelog

All notable changes to this project will be documented in this file.

## UNRELEASED

### Added

- **[added]**: new parameter (`min_reads_for_taxid`)to set a threshold for minimum number of reads assigned for a taxid to be considered
- **[added]**: new workflow (`SUBTYPE_AND_SEGMENT_FLU.nf`) attempts to retrieve the flu subtype and segment from kraken report file and populates the meta with these values
- **[added]**: new module (`retrieve_flu_subtype_and_segment.nf`) attempts to parse out the flu subtype and segment from the kraken report file and sets these values to Null if nothing retrieved

### Changed

- **[improvement]**: Taxid reference fasta files for consensus sequence are obtained from kraken database
- **[improvement]**: Channels now rely on Meta Mapping
- **[improvement]**: Output folder now have the following structure `output/<sample_id>/<taxid>`
- **[improvement]**: `write_manifest.py` relies on glob expression

### Removed

- **[Removed]**: writing manifest process removed from `SORT_READS_BY_REF.nf` 
- **[Removed]**: json resource files and fasta files provided on the repo

### Fixed

---
## [0.0.1] - 2023-12-01

This is the first prototype versioning. This pipeline 1) sort reads via Kraken and 2) generate consensus sequences using ivar.
