# Insertion mapper with BLAST

Short-reads generate DUP calls for genome duplications, which annotate the reference region that was duplicated.
Long-reads call duplications as large insertions (INS calls), which gives the location where the duplicated copy was
inserted and the sequence of the duplicated copy. Not all insertions can be represented as a duplicated reference locus,
but it is often necessary to determine the duplicated locus from long-read insertion calls, and mapping insertion
sequences with modern aligners often produces incomplete duplicated loci. For tandem duplications, the true breakpoint
is often shifted and embedded within the insertion sequence, and when mapping back to the reference, one side of the
duplication is lost and left unmapped.

To find duplication loci, this pipeline uses BLAST to seed alignments and several heuristics to locate a putative
duplicated locus for insertions. Better solutions will come from improved alignment and variant discovery methods,
however, this pipeline is intended to be used in the interim to improve callsets and aid method development.

## Requirements

### Docker and Singularity

A containerized version of this pipeline is not yet available, but is intended. All required tools and libraries will
be part of the embedded in the container.

### Base pipeline

Python 3+ with libraries:
* biopython
* numpy
* pandas
* pybedtools
* pysam
* snakemake

Other tools:
* NCBI BLAST
* BedTools
