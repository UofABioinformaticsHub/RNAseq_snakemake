# University of Adelaide Snakemake RNA-seq Analysis Pipeline

Snakemake is a workflow management system which can be used to create reproducible analysis pipelines. Here, we have developed a standardised pipeline for RNA-seq analysis, that takes paired raw sequence data as input and outputs count matrices from a range of quantification tools.

## Pipeline setup

The instructions below will help prepare your Phoenix environment to be able to run the pipeline, format your data and generate the few assistance files required by the pipeline. 

## Installing the Pipeline

To install the pipeline scripts, simply clone this repository:

```
git clone https://github.com/UofABioinformaticsHub/RNAseq_snakemake
```

## Setting up a Python Virtual Environment

#### Background on Python Virtual Environment

Python applications often require packages and modules (utilities/tools) that don't come with the standard python library. Often, software we want to use will require a specific version of python, python library or an accompanying piece of software. The reason for these dependencies can range from how the software was built to it requiring obsolete dependencies. 

Consequently, it's not always possible for the system wide python installation to meet the needs of every piece of software we want to use. A solution to this problem is to create and use a python virtual environment (PVE). PVEs are self contained directories that can be built with any python version, where we can install any software we desire. 
The benefit of PVEs is that they are independent of the parent system, meaning they do not intefere with system wide python packages or settings. It also means each python tool we want to use can have its own PVE, which will avoid any incompatibility issues that might otherwise be faced. 

#### Setting up a Python Virtual Environment

The Phoenix team have written a great [tutorial](https://wiki.adelaide.edu.au/hpc/index.php/Python_virtual_environment) on activating and using PVEs, which we recommend you follow. However, a breif summary of the set up process is below.

First, load the module for python 3.6.1 in a Phoenix session

```
$ module load Python/3.6.1-foss-2016b
```

Then, create the PVE in your `Fast` directory (or any directory you like really...), giving the PVE a sensible name that has meaning to its function. Snakemake requires python 3, so we will need to specify this when creating the environment. 

```
$ virtualenv --system-site-packages -p python3 $FASTDIR/path/to/desired/location/py_snakemake
```

That's it! The above command will create the directory `py_snakemake` which will have everything it needs to function as a PVE.

## Sample sheet 

The RNAseq pipeline uses what's called a `YAML` configuration file to provide information to the `snakefile`, the script that is executed. Configuration files are beneficial as they can be edited on the fly, with the resulting changes flowing on to the executable script. 
Here, we've written a bash script that generates a YAML configuration file when provided with a `.tsv` file containing sample information.
In the GitHub repository, there is a template samplesheet file called `samplesheet.tsv`. The column headers are described in greater detail below. The `RG_` prefix stands for `readgroup`, and is information that will be used to generate readgroups for the alignment files (BAM files).

- **FILENAME**
    - Full filenames of the fastq files, not including file extensions (e.g. DO NOT include \_R1.fastq.gz or \_R2.fastq.gz).
- **RG_SM** 
    - Shortened unique sample identifier which is used for BAM readgroups. This can often be a unique fragment of the filename above.
- **RG_ID** 
    - This column is to specify which group reads belong to. It is built from the header line of the fastq files.
- **RG_PU**
    - The platform unit is a combination of the flowcell barcode, the lane number and sample barcode.
- **RG_PL** 
    - Platform the sequences were generated on.
- **RG_LB** 
    - Library preparation method.
- **R1_extension**
    - This is the extension for the forward reads (_e.g. \_R1.fastq.gz, \_R1.fastq, \_R1.fq.gz etc..._).
- **R2_extension**  
    - This is the extension for the reverse reads .
- **PATH**
    - The file path to where the reads are located on Phoenix. Remember to include a final backslash (/) at the end of the path!

Much of how to generate the IDs for the readgroups is explained in detail at the [following link](https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups).
An example `sampleSheet.tsv` with three paired-end samples is presented below.

**FILENAME**|**RG\_SM**|**RG\_ID**|**RG\_PU**|**RG\_PL**|**RG\_LB**|**R1\_extension**|**R2\_extension**|**PATH**

:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
sample1\_filename|N9001|NB501008.4|NB501008HYV7VBGXY.4|ILLUMINA|made\_up\_LB|\_R1.fastq.gz|\_R2.fastq.gz|/fast/users/a1234567/data/fastq/
sample2\_filename|AX834|NB501008.4|NB501008HYV7VBGXY.4|ILLUMINA|made\_up\_LB|\_R1.fastq.gz|\_R2.fastq.gz|/fast/users/a1234567/data/fastq/
sample3\_filename|LY610|NB501008.4|NB501008HYV7VBGXY.4|ILLUMINA|made\_up\_LB|\_R1.fastq.gz|\_R2.fastq.gz|/fast/users/a1234567/data/fastq/

## Building the `YAMl` configuration file

Once the samplesheet has been filled out, building the configuration file is simply executing the `tsv2yaml.sh` script in the same directory as `sampleSheet.tsv`. 

```
$ ./tsv2yaml.sh sampleSheet.tsv
```

This will make the configuration file which will appear in the same directory as where `tsv2yaml.sh` was executed. 
The next step is to edit the file to include your specific file paths/data files you'll need. An example of the configuration file is shown below.

```
THREADS: 8
BASE: '/fast/users/a1234567/snakemake/full_pipeline'
RAW_FILE_PATH: '/fast/users/a1234567/snakemake/fastq'
FASTQC: '/0_FastQC'
QC_raw: '/FastQC_raw'
QC_trim: '/FastQC_trim'
...
TRAN: '/data/biohub/Refs/zebrafish/Danio_rerio.GRCz10.cdna.all.fa'
GTF: '/data/biohub/Refs/zebrafish/Danio_rerio.GRCz10.88.gtf'
...
SampleList:
  Sample_1_A:
    Sample: 'Sample_1_A'
    RG_SM: 'smp1_A'
    RG_ID: 'NB501008.4'
    RG_PU: 'NB501008HYV7VBGXY.4'
    RG_PL: 'ILLUMINA'
    RG_LB: 'made_up_LB'
    R1: '/fast/users/a1234567/snakemake/fastq/Sample_1_A_R1.fastq.gz'
    R2: '/fast/users/a1234567/snakemake/fastq/Sample_1_A_R2.fastq.gz'
```

When editing the configuration file, it is important not to change any of the indentation or syntax as this will affect how the information is read, likely resulting in a keyerror. 

## Editing the Phoenix batch script

There is a template Phoenix batch script in the GitHub repository called `phoenix_snakemake.sh`. This is like any other batch script that you use to submit a job to phoenix. A breakdown of the script is below.

#### Job parameters

```
#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=2:00:00
#SBATCH --mem=64GB
#SBATCH -o /fast/users/a1645424/snakemake/test_output/slurm/full_pipeline_%j.out
#SBATCH -e /fast/users/a1645424/snakemake/test_output/slurm/full_pipeline_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au
```

The first section is the job parameters, which tell Phoenix how many resources your job requires. If you are familiar with submitting jobs to Phoenix, you can edit this section to look however you like.
It's important that the number of threads specified in the batch script matches or exceeds those specified in the configuration file. If fewer threads are specified in the batch script, snakemake will account for this by using fewer threads, however the job will then be less efficient.
It's also important to allocate memory based on the most memory intensive task in the pipeline (likely the STAR indexing phase).
Hardset memory and CPU limits are a little rigid, as not all jobs are going to need the same resources. A SLURM configuration file is is in the works, which will be capable of allocating memory and CPU requirements based on the specific rule.

#### Modules 

```
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
```

Currently, the snakemake pipeline utilises the module system on Phoenix. The order of the modules is important, as loading multiple modules can result in software conflicts that affect the operation of loaded tools. 

#### The snakemake executable 

```
source /fast/users/a1645424/snakemake/py_snakemake/bin/activate
snakemake -R --cores 8 -s /fast/users/a1645424/snakemake/scripts/snakefile.py
deactivate
```

The source command activates the PVE that we generated above. Ensure that the path to the `activate` executable is correct, otherwise it will throw an error message. 
The snakemake command is next. Here, users can take liberty in how they want to execute the pipeline. A basic command has been included, which will simply execute the whole pipeline. Be sure to change the number of cores to a value that matches your config/allocated resources, and check the file path to the snakefile is correct (see note on cores/threads below). 
If users want to only run a specific rule, or perform a dry run (i.e. check that all the dependencies between rules are valid without actually executing the code), then it is simply a matter of editing this snakemake command. Check out the [snakemake documentation](http://snakemake.readthedocs.io/en/latest/executable.html) for executable options.
Finally, the `deactivate` command is closing the PVE after the snakemake pipeline has been run.

**NOTE:** 
The number of `--cores` and `--threads` mean sligtly different things. The number of `--cores` represents the total number of physical cores accessible by the pipeline. The number of `--threads` represents the number of threads dedicated to each job in the pipeline. 
For example, if you have four samples and want to run them all in parallel for a specific job (ensuring you have the memory overhead to do so), you could specify `--cores 16` and `--threads 4`. 
The main point is that the number of cores is the ultimate limiting factor. Snakemake will not run two samples in parallel for a job if the parallel thread count exceeds the specified cores count. 

## Submitting the pipeline to Phoenix

Once you've edited the `phoneix_snakemake.sh` script, you can submit the job to Phoenix using the command below 

`$ sbatch phoenix_snakemake.sh`

This will submit the job, giving you a job ID. You can then check the progress of the job using the following command:

`$ squeue -u a1234567`

If the pipeline fails, you can check what went wrong by checking the standard out and/or standard error files that we generated using the...

```
#SBATCH -o /.../.../full_pipeline_%j.out
#SBATCH -e /.../.../full_pipeline_%j.err
``` 

commands at the beginning of the batch script.

## Running specific rules

If you only want to run a specific rule/set of rules, there are a couple of ways to achieve this. The method I find most easy is to edit the `rule all` in the `snakefile.py`. The rule all dictates what actually needs to be run to generate the specified output files. If you only wanted to run FastQC, then you could specify only FastQC at the rule all like so:

```
QC_raw = config["BASE"] + config["FASTQC"] + config["QC_raw"]
QC_trim = config["BASE"] + config["FASTQC"] + config["QC_trim"]
rule all:
    input:
        QC_raw, QC_trim
```

In this example, the output of FastQC has been assinged to variables `QC_raw` and `QC_trim` to tidy the code. I find this method intuitive and simple as you are literally controlling what the pipeline generates. After editing the `snakefile.py`, execute the pipeline in the exact same was as shown above. 

<br>

A alternative method is to list the rules at the command line interface (i.e. edit the `snakemake -s /.../snakefile.py...etc` command in the `phoenix_snakemake.sh` script) that you want to re-run or only run. 
This could look something like the following if you wanted to re-run a specific rule...

```
snakemake -s snakemake.py -R rule_name_here
```

## Improvements/In development

- Report generation after FastQC and after sequence alignment
    - The code is there to generate the reports, the issue lies with rpy2 and the version of R it was compiled with.
- SLURM configuration file 
    - Simple YAML file that has a default resource configuration, as well as job specific resource configurations.
    - This will submit all samples for a specific snakemake job to phoenix in parallel (i.e. all samples will have their own Phoenix instance, resources etc...)

## Acknowledgements 

- The University of Adelaide Bioinformatics Hub: [GitHub Repo](https://github.com/UofABioinformaticsHub)
- Alastair Ludington: alastair.ludington@adelaide.edu.au
- Jimmy Breen: jimmy.breen@adelaide.edu.au
- Stephen Pederson: stephen.pederson@adelaide.edu.au
