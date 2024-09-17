# Software Requirements

Here is a list of all the softwares versions currently in use in the pipeline as set on each container.

- **Base container**

  - Samtools = `1.8`
  - BWA = `0.7.17`
  - biopython = `1.79`
  - pysam = `0.22.0`
  - pandas = `1.1.5`
  - matplotlib = `3.3.4`

- **Ivar container**
  - samtools = `1.11`
  - BWA = `0.7.17`

- **Kraken2Ref cotainer**:
  - pytest = `6.2.2`
  - importlib-resources = `5.1.0`
  - flake8 = `7.0.0`
  - pandas = `2.1.4`
  - cached-property = `1.5.2`
  - scipy = `1.12.0`
  - kraken2ref = `v2.0.0`

- **Kraken containers**
  - kraken2 = `v2.1.3`
  - kraken_tools = `v1.2`
  - biopython = `v1.81`

If not running under the provided containers, is strongly recommended to use the same versions described here.