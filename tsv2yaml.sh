#!/bin/bash -l

# Jimmy Breen & Alastair Ludington
# 170228

# This script converts the Biohub RNAseq samplesheet to a yaml format
# Samplesheet looks like this:

# Usage:
# pipeline_tsv2yaml.sh [SampleSheet.tsv]

# Read in file from commandline argument
sed 1d $1 > tmp.tsv

# replace all instances of '...' with a filepath/directory structure
# Select the BAM key corresponding to the aligner to be used
# Change "config...(dot)yaml" to whatever you want to call it
cat > config.yaml <<__SCRIPT__
## Directory paths
BASE: '/path/to/desired/output/destination'
QC_raw: '/0_FastQC_raw'
QC_trim: '/0_FastQC_trim'
TRIM: '/1_AdapterRemoval'
# BAMS: '/2_Hisat2_sambamba_merged'
# BAMS: '/2_STAR_sambamba'
unmapped: '/unmappedReads'
BAMS_RG: '/2_bams_read_group'
QUANT_SAL: '/3_salmon'
QUANT_FC: '/3_featureCounts'
ASS: '/4_stringtie'
KAL: '/5_kallisto'
LOG: '/logs'

## Data objects (afw)
TRAN: '.../cdna.fa.gz'
INDEX_hst2: '...'
GTF: '...'
KAL_idx: '.../.idx'

## Function pathways (afw)
AdapterRemoval: '.../AdapterRemoval'
featurecounts: '.../featureCounts'
stringtie: '.../stringtie'
salmon: '.../salmon'
hisat2: '.../hisat2'
sambamba: '/localscratch/Programs/sambamba'
picard: '.../picard.jar'
kallisto: '.../kallisto'

## Adapter removal
basename: '...'
minquality_AR: 10
minlength_AR: 25
AR_threads: 1
DISCARD: '/AdapterRemoval/excessReads'

## Fastqc - parameters apply to both raw and trimmed fastq files
format: 'fastq'
kmer: 9
FQC_threads: 1

## Salmon quantifiation
     # INDEX_DIR - change the file path to where the index will be built (e.g. in QUANT_SAL)
INDEX_DIR: '.../Index'
bootstrap: 100
SAL_threads: 1

## hisat2_sambamba
quality: 30
H2_threads: 1

## STAR aligner
genomeDir: '/data/biohub/Refs/zebrafish/STAR_index'
compress: 5
ST_threads: 1

## FeatureCounts
minquality_FC: 10
strandness: 1
FC_threads: 1

## Stringtie
minlength_ST: 30
ST_threads: 1

## kallisto
bootstrap: 100
KAL_threads: 1

## Edit this section
SampleList:
__SCRIPT__

while read line
 do

################ ATTENTION ################
# If you have changed the column names/number of columns in the sample sheet (.tsv),
# you will need to alter the variable selection and dictionary structure below. 


     # Cutting samplesheet columns and assigning them to variables
     Sam=$(echo "$line" | cut -d $'\t' -f1)
     Plat=$(echo "$line" | cut -d $'\t' -f2)
     Rep=$(echo "$line" | cut -d $'\t' -f3)
     Treat=$(echo "$line" | cut -d $'\t' -f4)
     Idx=$(echo "$line" | cut -d $'\t' -f5)
     Suff=$(echo "$line" | cut -d $'\t' -f6)
     Path1=$(echo "$line" | cut -d $'\t' -f7)
     Path2=$(echo "$line" | cut -d $'\t' -f8)

     # Creating dictionary for each sample - from samplesheet
     echo "  "$Sam":" >> config.yaml
     echo "    Sample: "\'$Sam\'"" >> config.yaml
     echo "    Platform: "\'$Plat\'"" >> config.yaml
     echo "    Replicate: "\'$Rep\'"" >> config.yaml
     echo "    Treatment: "\'$Treat\'"" >> config.yaml
     echo "    Index: "\'$Idx\'"" >> config.yaml
     echo "    Suffix: "\'$Suff\'"" >> config.yaml
     echo "    R1: "\'$Path1\'"" >> config.yaml
     echo "    R2: "\'$Path2\'"" >> config.yaml
# Now we have to input the samplesheet we want to parse but first remove first line
done < tmp.tsv

rm tmp.tsv
