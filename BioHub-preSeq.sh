#!/bin/bash

module load AdapterRemoval/2.2.1-foss-2016b
module load Subread/1.5.2-foss-2016b 		# module load featurecounts WEHI
module load StringTie/1.3.3-foss-2016b
module load Salmon/0.8.2
module load HISAT2/2.0.5-foss-2016b
module load sambamba/0.6.6-foss-2016b
module load kallisto/0.43.1-foss-2016b
module load snakemake/3.12.0-foss-2016b-Python-3.6.1

module load picard/2.2.4-Java-1.8.0_121
java -jar /apps/software/picard/2.2.4-Java-1.8.0_71/picard.jar
