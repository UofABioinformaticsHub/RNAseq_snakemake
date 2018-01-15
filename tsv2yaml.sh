#!/bin/bash -l

# Jimmy Breen & Alastair Ludington
# 170228

# This script converts the Biohub RNAseq samplesheet to a yaml format
# Samplesheet looks like this:

# Usage:
# pipeline_tsv2yaml.sh [SampleSheet.tsv]

# Read in file from commandline argument
sed 1d $1 > tmp.tsv

## << denotes a here document - 
cat > config.yaml <<__SCRIPT__ 
THREADS: 8

## Directory paths
BASE: '/path/to/parent/output/directory'
RAW_FILE_PATH: '/path/to/raw/fastq/files'

FASTQC: '/0_FastQC'
QC_raw: '/FastQC_raw'
QC_trim: '/FastQC_trim'

AR: '/1_AdapterRemoval'

ALIGN: '/2_STAR'
BAMS: '/bams'
RG_PICARD: '/bams_readgroups'

QUANT_SAL: '/3_salmon'

QUANT_FC: '/4_featureCounts'

QUANT_STR: '/5_stringtie'

QUANT_KAL: '/6_kallisto'

LOG: '/logs'

## Data objects (afw)
TRAN: '/path/to/cdna/transcripts/fasta.fa'
GTF: '/path/to/organism/gtf/file.gtf'

## Adapter removal
minquality_AR: 10
minlength_AR: 25
AR_threads: 1
DISCARD: '/discarded_reads' ## not currently in use

## Fastqc - parameters apply to both raw and trimmed fastq files
format: 'fastq'
kmer: 9

## STAR index
overhang: 149

## STAR mapping
BAMcompression: 1 ## Value between 1-10 where 1 is least compressed
qualityFilter: 30

## Salmon quantifiation
bootstrap: 100

## FeatureCounts
minquality_FC: 10
strandness: 1

## Stringtie
minlength_ST: 30

## kallisto
bootstrap: 100

SampleList:
__SCRIPT__

while read line; do

     # Parsing samplesheet columns
     SAM=$(echo "$line" | cut -d $'\t' -f1)
     RG_SM=$(echo "$line" | cut -d $'\t' -f2)
     RG_ID=$(echo "$line" | cut -d $'\t' -f3)
     RG_PU=$(echo "$line" | cut -d $'\t' -f4)
     RG_PL=$(echo "$line" | cut -d $'\t' -f5)
     RG_LB=$(echo "$line" | cut -d $'\t' -f6)
     R1_extension=$(echo "$line" | cut -d $'\t' -f7)
     R2_extension=$(echo "$line" | cut -d $'\t' -f8)
     Path=$(echo "$line" | cut -d $'\t' -f9)

     ## Creating YAML keys using parsed samplesheet variables
     echo "  "$SAM":" >> config.yaml
     echo "    Sample: "\'$SAM\'"" >> config.yaml
     echo "    RG_SM: "\'$RG_SM\'"" >> config.yaml
     echo "    RG_ID: "\'$RG_ID\'"" >> config.yaml
     echo "    RG_PU: "\'$RG_PU\'"" >> config.yaml
     echo "    RG_PL: "\'$RG_PL\'"" >> config.yaml
     echo "    RG_LB: "\'$RG_LB\'"" >> config.yaml
     echo "    R1: "\'$Path$SAM$R1_extension\'"" >> config.yaml
     echo "    R2: "\'$Path$SAM$R2_extension\'"" >> config.yaml
     
# Now we have to input the samplesheet we want to parse but first remove first line
done < tmp.tsv

rm tmp.tsv
