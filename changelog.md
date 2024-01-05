# Changelog

All notable changes to this project will be documented in this file.

## UNRELEASED

### Changed

- **[improvement]**: Channels now rely on Meta Mapping
- **[improvement]**: Output folder now have the following structure `output/<sample_id>/<taxid>`
- **[improvement]**: `write_manifest.py` relies on glob expression

### Removed

- **[Removed]**: writing manifest process removed from `SORT_READS_BY_REF.nf` 

### Added

### Fixed

---
## [0.0.1] - 2023-12-01

This is the first prototype versioning. This pipeline 1) sort reads via Kraken and 2) generate consensus sequences using ivar.