# Insertion mapper with BLAST

This pipeline is under active development as of 2023-12-13. New features are being added, it is being tested, and
documentation is being completed. A more mature version should be released within a few weeks.

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


## Pipeline structure

This documentation will refer to two specific directories:
* PIPELINE_DIR: Directory where the pipeline files are found. It will contain `libdupmap`, `rules`, and `Snakefile`.
* RUN_DIR: Directory where the pipeline is run.

Generally, PIPELINE_DIR and RUN_DIR should not be the same location if the pipeline is used for several projects.

The pipeline takes two input files:
* Sample table: Points to input files for each sample.
* Configuration file: Configuration parameters and pipeline tuning.

These inputs are described below.

## Docker and Singularity

A containerized version of this pipeline is not yet available, but is intended. All required tools and libraries will
be part of the embedded in the container. The full PIPELINE_DIR and all dependencies will be embedded in the container
once it is released.

## Native install

Until the Docker and Singularity packages are released, the package must be run from a native install. This requires
pulling the DupMapper code from GitHub and installing all dependencies.

### Clone

Change to the directory where the code will be installed. The directory `git clone` creates will become the
PIPELINE_DIR.

git clone --recursive https://github.com/BeckLaboratory/dupmapper.git

### Dependencies

Python 3+ with libraries:
* biopython
* numpy
* pandas
* pybedtools
* pysam
* snakemake
* intervaltree

Other tools:
* NCBI BLAST
* BedTools
* bcftools (for VCF input)


## Sample table

The sample configuration file tells the pipeline where to find input files for each sample. By default, this table is
called `samples.tsv`. If `samples.tsv` is missing, it will search for `samples.xlsx`, and if present, the first sheet
from the Excel file is read. If the `sample_table` configuration option is set (see below), the file it points to is
read and the default files are ignored if present. 

Table must have two columns:
* SAMPLE: Sample name. Each row must have a unique name, no blank names are allowed.
* DATA: Path to input data.

Each entry may point to a VCF or BED file.

### VCF input

A VCF file should contain one or more SV insertions (50+ bp by default). VCF records can either have the full SV
sequence in the REF and ALT columns, or it can use symbolic ALTs (e.g. "<INS>") if the sequence is in the INFO field.
By default, the sequence is expected as a "SEQ" INFO attribute, but the name of the field in INFO is tunable (see
"vcf_seq" configuration option).

If the ID is missing, one will be assigned for each record.

The DATA column for VCF input should be a path to the input VCF file.

### BED input

Several pipelines represent variants as BED files with standard fields, which can be pulled directly from SV-Pop or PAV.
Input BED files must have a set of standard fields:

* #CHROM: Chromosome name.
* POS: Reference start position.
* END: Reference end position.
* ID: A unique ID for each variant.
* SVTYPE: SV type. Entries where SVTYPE is not "INS" will be dropped.
* SVLEN: Variant length. 
* SEQ: Optional, variant sequence. If missing, sequence is found in a FASTA file (see below).

POS and END are in BED 0-based half-open coordinates.

If the BED file does not contain the variant sequence for each record in the SEQ field, then an external FASTA file
must be provided containing the sequence for each variant. The ID for records in the BED file must match FASTA record
names. The FASTA filename must also appear in the DATA field with BED filename separated by
a semicolon (;). For example, DATA="/path/to/sample/sv_ins.bed.gz;/path/to/sample/fasta/sv_ins.fa.gz"). If both the
SEQ field is present and a FASTA file is specified, the FASTA file will be used to locate sequences.

If the variants are pulled from SV-Pop, the pipeline can find the location of the FASTA file corresponding to the BED
file automatically without needing to concatenate the FASTA filename with the BED filename in DATA. It does this by
searching for a file pattern. If the BED file is "/path/to/sample/bed/FILENAME.bed.gz", then it searches for
"/path/to/sample/bed/fa/FILENAME.fa.gz". If the SEQ field is present, this path is not searched.

## Configuration file

The configuration file is located in the RUN_DIR and is called `config.yaml` by default. The file format is YAML, and
this documentation will give examples for its format. If `config.yaml` is missing but `config.json` is present, it
will be read as the configuration file in JSON format.

A minimal configuration `config.yaml` file has "reference" and "blast_db" parameters set:

```
---

reference: /path/to/ref/hg38/no_alt/hg38.no_alt.fa.gz
blast_db: /path/to/ref/hg38/no_alt/blastdb/2.13.0/hg38_no_alt
```

Required parameters:
* reference [str]: Path to the reference FASTA file. May be gzipped.
* blast_db [str]: Prefix of the BLAST database, e.g. PREFIX + ".nsq" is the BLAST sequence file. PREFIX may be a full
  path to the database files (e.g. "/path/to/blastdb/hg38" where "/path/to/blastdb/hg38.nsq" is one of the database
  files).

BLAST parameters:
* blast_param [str, "-word_size 16 -perc_identity 95"]: BLAST parameters.
* blast_range [str, "50:30k"]: Attempt BLAST for variant sizes in this range (min:max).

Alignment parameters:
* align_range [str: "50:"]: Attempt alignment for variant sizes in this range (min:max).
* aligner [str, "minimap2"]: 
* align_params [str, "-x asm5"]

Reference masking parameters:
* mask_merge_dist [int, 20]: 
* mask_bed_list [list(str)]: List of one or more BED files containing additional masked loci.
* ref_soft_mask [bool, true]: Use reference soft-masking to filter BLAST hits (recommended).

Input parameters
* sample_table [str]: Sample table filename. If absent, a sample table is searched for (see "Sample configuration").
* vcf_filter_gt [bool, True]: Filter on the genotype column. If a VCF record does not have at least one alternate 
  allele for a record, it is discarded. Leave on for multi-sample VCFs.
* vcf_pass [list, ["PASS", "."]]: A list of strings of acceptable FILTER column values for records. If the value is
  a string, then only FILTERs with that value are accepted. If the parameter is blank, then all values are accepted.
* vcf_strict_sample [bool, False]: For single-sample VCFs, require the sample column name to match the sample being
  processed. Turning this on for single-sample VCFs where the sample name this pipeline sees should match the VCF
  sample name will enable a check for incorrect sample input.
* vcf_seq [str, "SEQ"]: For VCF files with symbolic ALTs (e.g. "<INS>" instead of full sequence in REF & ALT columns),
search for the SV sequence in this info field.


Additional pipeline control:
* shell_prefix [str]: Pre-pend shell commands (such as BLAST commands) with this string. May be used to prepare the
  shell, for example, loading software modules.
  * BASH strict mode "set -euo pipefail" is prepended if missing, and a semi-colon is added to the end of the string
    if missing.
* temp_dir [str, "temp/working_temp"]:
* partitions [int, 80]:
* td_pad [float, 0.1]: 
* td_pad_min [int, 100]: 


### Advanced configuration

On the command-line, a path to a specific configuration file can be given (e.g. `--config config=/path/to/config.yaml`).
The configuration file must end with ".yaml" or ".json" (not case sensitive) and will be read as a YAML or JSON file,
respectively.

## Running the pipeline

This section assumes that the configuration file and sample table (`config.yaml` and `samples.tsv` by default) are
already in the RUN_DIR. See instructions above to set those up.

### Docker and Singularity

Coming soon.

### Native install

To run the native install, link `Snakefile` and `profiles` (optional, see below) from PIPELINE_DIR to RUN_DIR.

```
cd RUN_DIR
ln -s PIPELINE_DIR/Snakefile ./
ln -s PIPELINE_DIR/profiles ./
```

Load any environment modules or conda environments needed to place the Python and tool dependencies in the current
environment (i.e. Python 3, bcftools, BEDTools, etc) and change to the RUN_DIR.


#### Run locally

To run locally on the current machine:

```
snakemake --profile profiles/local --cores 32
```

* `--profile profiles/local`: Run the default local profile. See "Profiles" below.
* `--cores 32`: Consume up to 32 cores.

The `-n` runs through the pipeline without executing anything. It can be helpful for checking that input files are
being read properly.


#### Run distributed

The pipeline contains a default Slurm profile.

```
snakemake --profile profiles/slurm -j 20
```

* `--profile profiles/slurm`: Run the default slurm profile. See "Profiles" below.
* `-j 20`: Run 20 jobs simultaneously.


### Profiles

Profiles setup resources for the pipeline and control distributing jobs over a cluster. The pipeline has two default
profiles, one for running on the current machine (`profiles/local`) and one for distributing jobs over a cluster
(`profiles/slurm`).

Snakemake profile documentation can be found here:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles

To make changes to the profile, copy instead of linking the profiles directory (`unlink profiles` to remove the link if
it was already created):

For example, to customize the slurm profile:
```
mkdir profiles
cp -r PIPELINE_DIR/profiles/slurm profiles/
```

Then edit `profiles/slurm/config.yaml` for your environment.
