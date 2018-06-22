## Overview
This repository contains benchmark scripts for transcript quantification tools. Here is the list of the tools that used in this benchmark. If you want to add your tool into this benchmark, please just submit an github issue. 

- Bowtie2 + RSEM (Bowtie2 version v2.3.4.1, RSEM version v1.3.0)
- STAR + RSEM (STAR version 2.6.0c, RSEM version v1.3.0)
- Kallisto (Version 0.44.0)
- Salmon (Version 0.9.1)
- Hera (Version 1.2)
- Nora (Version 1.0.3 https://nora.bioturing.com/) 

## Data sets: 

We use simulated benchmark data in Kallisto's paper (Bray et al., 2016). In particular, this data set ''contains 20 RNA-Seq simulations generated with the program RSEM. The transcript abundances and error profiles for the
simulated data were based on the quantification of sample NA12716_7 from the 
GEUAVDIS dataset, and to accord with GEUVADIS samples the simulations consisted
of 30 million reads.'' Below are the data preparation steps. 

### Create RSEM index 
Install RSEM
```
$cd 
$wget https://github.com/deweylab/RSEM/archive/v1.3.0.tar.gz
$tar -xvf v1.3.0.tar.gz
$cd RSEM-1.3.0
$make
$export PATH=$PATH:~/RSEM-1.3.0
```
Download the reference and gene anotation files. 
```
$cd
$mkdir simulation
$mkdir simulation/data
$cd simulation/data
$wget ftp://ftp.ensembl.org/pub/release-80/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
$wget ftp://ftp.ensembl.org/pub/release-80/gtf/homo_sapiens/Homo_sapiens.GRCh38.80.gtf.gz
$gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
$gunzip Homo_sapiens.GRCh38.80.gtf.gz
```
Create index files
```
$cd ..
$mkdir index
$rsem-prepare-reference --gtf data/Homo_sapiens.GRCh38.80.gtf data/Homo_sapiens.GRCh38.dna.primary_assembly.fa ./index/grch38
```


We use Kallisto's benchmark sing 20 simulated data sets generated from Kallisto paper (using this script), and the most recent benchmark data from SMC-RNA DREAM challenge.



tool calculates the correlation between two
transcript expression files. This tool will output some stats:
  - Spearman: spearman correlation between tpm values
  - Pearson: pearson correlation between tpm values
  - Log-pearson: pearson correlation between log-transformed tpm values
  - MAE(asinh): mean absolute error of asinh-transformed tpm value
  - False positive: the number of unexpressed transcripts but predicted to be expressed by the program
  - False negative: the number of expressed transcripts but predicted to be unexpressed by the program
  - Max false neg: the maximum tpm value of the transcripts but predicted to be unexpressed by the program
  - Max false pos: the maximum predicted tpm value of the unexpressed transcripts

Spearman, pearson, log-pearson is calculated base on SMC-RNA-Challange evaluation script (https://github.com/Sage-Bionetworks/SMC-RNA-Challenge)

## Install

```shell
git clone https://github.com/bioturing/nora-benchmark
cd nora-benchmark
g++ benchmark.c -O2 -o benchmark
```

## Usage

```shell
Usage: ./benchmark TRUTH_FILE INPUT_FILE
```

```shell
Truth file and input must be a tsv (tab-separated) file — with no header. Containing two
columns mapping of each transcript present in the reference to the corresponding tpm values
(the first column is a transcript and the second is the corresponding tpm value).

ENST00000000233    28.410000
ENST00000000412    0.000000
ENST00000000442    0.000000
...
```

```shell
You can use this command to convert Nora output file to this format:
    tail -n +2 nora_out.tsv | awk '{printf "%s\t%s\n", $1, $8}' > input.tsv
```

# Contacts

Please report any issues directly to the github issue tracker. Also, you can send your feedback to support@bioturing.com
