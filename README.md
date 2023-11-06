# Nanamp, another workflow to analyse Oxford Nanopore SSU amplicon sequencing data 

## Why Nanamp?

This workflow has been developed to easily assigned amplicon from Oxford Nanopore sequencing using IDTAXA. This approach is slower than other approaches commonly used to assigned long read amplicons such as minimap and kraken but return much less false positives.

## The workflow

### Overview

This workflow is composed of four main steps. Primers are first removed using cutadapt, sequence are then filtered based on their quality and length using VSEARCH. Filtered sequences are clustered using VSEARCH's command `--cluster_fast` and clusters consensus sequences are taxonomically assigned using IDTAXA. In the end and OTU table with taxonomy assignments is returned.

## Getting started

### Before running the workflow

To run this workflow you need to have Nextflow installed and either Docker or Singularity, that's it.

### Running the workflow

You need first to have your files demultiplexed, which may already be the case. Once it is done, you need to assemble a two columns table making the correspondence between sample names and file names. An example below:

```
sample_1 data/sample1.fastq.gz
sample_2 data/sample2.fastq.gz
```

Fields (columns) are separated by a single space. By default this table is in a text file named `manifest.txt` located at the root of your working directory. We will see later how to specify another path for the file containing this table.

An example of how to run the workflow with the default parameters using docker:

```
nextflow run nhenry50/nanamp -r main -profile docker
```

The default parameters may not be suited for your data. You can override default parameters by indicating new values in a yaml file. Below are listed the default parameters, suitable for 16S amplicons amplified with the 27F/1492R primer pair:

```
# Input options

manifest: 'manifest.txt'

# Output options

outdir: 'outputs'

# Primer trimming options

primer_f: 'AGAGTTTGATCMTGGCTCAG'
primer_r: 'AAGTCGTAACAAGGTAACC'
cutadapt_error: 0.2

# Filtering options

min_length: 1400
max_length: 1800
maxee_rate: 0.2

# Clustering options

clusteringid = 0.97
clusteringiddef = 0

# Taxonomic assignment options

refdb = 'silva'
fastachunks = 1000
idtaxa_thresh = 50

# Max resource options

max_memory = '128.GB'
max_cpus = 16
max_time = '240.h'
```

You can simply paste these parameters in a file named `params.yaml` and change the ones you want. After that, you can run the workflow with your parameters with the option `-params-file`:

```
nextflow run nhenry50/nanamp -r main -params-file params.yaml
```

## Parameters

### Input options

`manifest`: file with the list of fastq files to analyse and their corresponding file names. Default: `manifest.txt`

### Output options

`outdir`: directory where to save the analyses results. Sub-directories will be created for each step of the analysis producing results. Default: `outputs`

### Primer trimming options

`primer_f`: forward primer sequence. Default: `AGAGTTTGATCMTGGCTCAG`

`primer_r`: reverse complement of the reverse primer sequence. Default: `AAGTCGTAACAAGGTAACC`

`cutadapt_error`: maximum error rate allowed to find primers using cutadapt (-e option). See cutadapt documentation for more information (https://cutadapt.readthedocs.io/en/stable/index.html). Default: `0.2`

### Filtering options

`min_length`: sequence minimum length. Default: `1400`

`max_length`: sequence maximum length. Default: `1800`

`maxee_rate`: maximum average expected error. See vsearch documentation (FASTA/FASTQ/SFF file processing options) for more details. Default: `0.2`.

### Clustering options

`clusteringiddef`: threshold for clustering value ranging from 0 to 1. See vsearch documentation for more details. Default value: `0.97`

`clusteringiddef`: type of distance measure used for clustering. Default value: `0`. From vsearch documentation:

```
Change the pairwise identity definition used in --id. Values accepted are:

      0.  CD-HIT definition: (matching columns) / (shortest sequence length).

      1.  edit distance: (matching columns) / (alignment length).

      2.  edit distance excluding terminal gaps (default definition for --id).

      3.  Marine  Biological  Lab  definition counting each gap opening (internal or terminal) as a single mis‐
          match, whether or not the gap was extended: 1.0 -  [(mismatches  +  gap  openings)/(longest  sequence
          length)]

      4.  BLAST definition, equivalent to --iddef 1 for global pairwise alignments.

The option --userfields accepts the fields id0 to id4, in addition to the field id, to report the pairwise iden‐
tity values corresponding to the different definitions.
```

### Taxonomic assignment options

`refdb`: reference database used for taxonomic assignment. Two choices, `silva` or `pr2`

`fastachunks`: number of sequences to assign per process. No need to change to change this parameter.

`idtaxa_thresh`: numeric specifying the confidence at which to truncate the output taxonomic classifications. Default value: `50`, same as the default value of the function `IdTaxa()` from the R package `DECIPHER`. more information can be found in the IdTaxa documentation.

### Max resource options

These parameters are the same as the one used in nf-core workflows:

`max_memory`: Maximum number of CPUs that can be requested for any single job.

`max_cpus`: Maximum amount of memory that can be requested for any single job.

`max_time`: Maximum amount of time that can be requested for any single job.
