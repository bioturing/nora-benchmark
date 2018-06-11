## Overview

Developed by BioTuring (www.bioturing.com), this tool calculates the correlation between two
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
