#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=2:00:00
#SBATCH --mem=64GB
#SBATCH -o /path/slurm/output/directory/slurm_%j.out
#SBATCH -e /path/slurm/output/directory/slurm_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=user.email@adelaide.edu.au


## Modules needed for the pipeline to run
## The order of these is important as certain modules will replace compilers/software versions. 
## Kallisto is currently loaded here and in the shell script in the snakefile. Probably can remove it from here. 

module load Python/3.6.1-foss-2016b
module load AdapterRemoval/2.2.1-foss-2016b
module load fastqc/0.11.4
module load Subread/1.5.2-foss-2016b
module load StringTie/1.3.3-foss-2017a
module load Salmon/0.8.2
module load STAR/2.5.3a-foss-2016b
module load sambamba/0.6.6-foss-2016b
module load picard/2.2.4-Java-1.8.0_121
module load kallisto/0.43.1-foss-2017a
module load SAMtools/1.3.1-foss-2016b

## Activating python virtual environment
source /path/to/python/virtal/environment/bin/activate 

## Snakemake command
## Can alter this to do whatever you please
## The `-R` command is included to prevent the 'snakefile already writing to the directory...blah blah blah' error that can occur if you run the same snakefile in two different shell sessions.
	## I.e. if you re-run this script due to an error. 
snakemake -R --cores 16 -s /path/to/snakefile.py

## deactivating python virtual environment
deactivate 
