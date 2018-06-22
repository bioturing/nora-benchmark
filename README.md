## Overview
This repository contains benchmark scripts for transcript quantification tools. If you want to add your tool or any additional data sets into this benchmark, please create a pull request or submit a github issue. 

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

### Generate Simulated Data
# Install RSEM
```shell
cd 
wget https://github.com/deweylab/RSEM/archive/v1.3.0.tar.gz
tar -xvf v1.3.0.tar.gz
cd RSEM-1.3.0
make
export PATH=$PATH:~/RSEM-1.3.0
```
Download the reference and gene anotation files. 
```shell
cd
mkdir simulation
mkdir simulation/data
cd simulation/data
wget ftp://ftp.ensembl.org/pub/release-80/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/gtf/homo_sapiens/Homo_sapiens.GRCh38.80.gtf.gz
gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.80.gtf.gz
```
Create index files
```shell
cd ..
mkdir index
rsem-prepare-reference --gtf data/Homo_sapiens.GRCh38.80.gtf data/Homo_sapiens.GRCh38.dna.primary_assembly.fa ./index/grch38
```

Download the transcript abundances and error profiles of the simulated data from Kallisto's github. These were based on the quantification of sample NA12716_7 from the GEUAVDIS dataset.
```shell
cd ..
mkdir NA12716_7
cd NA12716_7
wget https://raw.githubusercontent.com/pachterlab/kallisto_paper_analysis/nbt/simulations/NA12716_7/rsem/out.stat/out.model
wget https://raw.githubusercontent.com/pachterlab/kallisto_paper_analysis/nbt/simulations/NA12716_7/rsem/out.isoforms.results
```

Generate simulated data 

```shell
cd ..
for NUM in {1..20};
do
        rsem-simulate-reads index/grch38 NA12716_7/out.model NA12716_7/out.isoforms.results 0.0 30000000 sim_$NUM --seed $NUM
done
```
Download, compile, and put the correct binary versions of Salmon, Kallisto, Hera, Nora, Bowtie2, STAR into the correct executable paths. 

Generate index files for each tool:

```shell
GENOME_GTF=data/Homo_sapiens.GRCh38.80.gtf
GENOME_FASTA=data/Homo_sapiens.GRCh38.dna.primary_assembly.fa

mkdir ./index/star_rsem
rsem-prepare-reference --gtf $GENOME_GTF --star -p 32 $GENOME_FASTA ./index/star_rsem/grch38
mkdir ./index/nora
Nora index -g $GENOME_FASTA -t $GENOME_GTF -p genome -o ./index/nora/
mkdir ./index/hera
hera_build --fasta $GENOME_FASTA --gtf $GENOME_GTF --outdir ./index/hera/grch38
mkdir ./index/bowtie2_rsem
rsem-prepare-reference --gtf $GENOME_GTF --bowtie2 -p 32 $GENOME_FASTA ./index/bowtie2_rsem/grch38
salmon index --index ./index/salmon --transcripts ./index/rsem/grch38.transcripts.fa
kallisto index -i ./index/kallisto ./index/rsem/grch38.transcripts.fa 
```

Run paired-end mode
```shell
for i in {1..20}
do
	if [ -f sim_${i}_1.fq ] && [ -f sim_${i}_2.fq ] ; then
		./runPairedEnd.sh sim_${i}_1.fq sim_${i}_2.fq ./result/sim${i}
	fi
done
```
Run single-end mode
```shell
for i in {1..20}
do
	if [ -f sim_${i}_1.fq ] ; then
		./runSingleEnd.sh sim_${i}_1.fq ./result/sim_single${i} 180.5069 52.54896
done
# 180.5069 and 52.54896 are the mean fragment length and standard deviation, respectively. 
````
Given the known truth tpm, we use the following metrics for evaluations:

  - Spearman: spearman's rank correlation between tpm values
  - Pearson: pearson correlation between tpm values
  - Log-pearson: pearson correlation between log-transformed tpm values
  - MAE(asinh): mean absolute error of asinh-transformed tpm value
  - False positive: the number of unexpressed transcripts but predicted to be expressed by the program
  - False negative: the number of expressed transcripts but predicted to be unexpressed by the program
  - Max false neg: the maximum tpm value of the transcripts but predicted to be unexpressed by the program
  - Max false pos: the maximum predicted tpm value of the unexpressed transcripts
  
 To evaluate, we first have to prepare the truth tpm and estimated tpm files. Each file has 2 columns, transcript id and tpm value without headers.

Compile benchmark.cpp then run ./benchmark <truth_file> <estimated_file>

Spearman, pearson, log-pearson is calculated base on SMC-RNA-Challange evaluation script (https://github.com/Sage-Bionetworks/SMC-RNA-Challenge)

## Install

```shell
git clone https://github.com/bioturing/nora-benchmark
cd nora-benchmark
g++ benchmark.cpp -O2 -o benchmark
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
