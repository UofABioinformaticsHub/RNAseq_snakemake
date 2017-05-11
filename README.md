# __RNAseq_snakemake__

This is an RNAseq analysis pipeline written in snakemake for the UofA Bioinformatics Hub. An outline of the pipeline is presented below

- PART0: Targets for pipeline
- PART1: AdapterRemoval fastq
- PART2: QC raw fastq
- PART3: QC trimmed fastq
- PART4: Hisat2/Sambamba
- PART4.5: STAR/Sambamba
- PART5: Generating readgroups
- PART6: Quantification Salmon
- PART7: FeatureCounts
- PART8: Stringtie assembly part 1
- PART8.5 Stringtie making part 2
- PART9: Kallisto

### __Installing snakemake__

Snakemake is available via module on phoenix (need to check specifics). To install Snakemake on a personal computer, see the [installation page](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) from the snakemake documentation.

### __Usage__

Clone the repository to obtain the `snakefile`, `sampleSheet` and configuration file compiler `tsv2yaml`. Fill in the sample sheet with sample information, being sure to keep the tab separated format. Columns can be added/removed from the sample sheet (or an Illumina samplesheet can be used in replacement in tsv format) but will result in mandatory edits to the the config compiler and snakefile. Specify directory structure and program parameters in the `tsv2yaml.sh` script.

To generate the configuration file, run the `tsv2yaml.sh` script as below:

```./tsv2yaml.sh sampleSheet.tsv```

The default output will be a configuration file in yaml format labelled `config.yaml`. Any changes made to key values (any string before a colon '`:`''), dictionary structure (tab indentation for keys) or config file name are likely to result in incompatibilities between the snakefile and the config file. Be sure to go through the snakefile and make the necessary amendments to the affected fields (e.g. config[""} in the snakefile).


##### _Standard execution_

To run the snakefile, specify the snakemake environment and provide the snakefile as below:

```snakemake -s snakefile```

The `-s` argument is used to specify a snakefile. If the snakefile is in the directory where the snakemake executable is run and does not have any extensions, the `-s` argument is not necessary.


##### _Dry run_

A dry run tests whether the rules in the snakefile have been defined correctly, without generating any output. Dry runs can be executed with the `-n` argument:

```snakemake -s snakefile -n

If there are errors in the snakefile, whether it's a typo in the code, missing input files or undefined variables, a dry run will provide a detailed message of what went wrong and on what line.

##### _Overriding thread usage and parallelisation_

From the commandline, the number of jobs to be run in parallel can be set using the `-j` argument:

```snakemake -j 5 -s snakefile```

The above command will execute five instances of the snakefile in parallel. It is important to remember how many threads have been specified, as jobs $\times$ threads == total threads used. For example, running 5 jobs in parallel with 10 threads each is 50 threads in total.

It is also possible to adjust the number of cores used for each job using the `--core` argument:

```snakemake --core 10 -s snakefile```

By specifying `--core`, thread values set throughout the snakefile (either hard coded or specified from a config file) will now be set to 10 threads. This is a global alteration and will be applied to all jobs!

##### _Running a specific rule_

Snakemake will identify whether a rule needs to be re-run, such as when the input files are newer than the output files or a new output file is ammended into in the snakefile. However, just altering rule parameters is not enough to trigger a re-run of a rule. This is because snakemake only looks to see if the target files are there, not which parameters those files were made under.

To manually re-run a rule, use the `-R` argument followed by the rule name:

```snakemake -s snakemake -R madeUpRule```

The code above will force the re-execution of the rule `madeUpRule`, generating new output.

<br>

There are many more snakemake executables which you can read about [here](https://snakemake.readthedocs.io/en/latest/executable.html) or by running `snakemake -h` at the command line.

### __Contributing__

If you want to edit your own version of the snakemake pipeline, fork your own branch and go ham! The master branch will be updated with improved software as necessary, however a stable (or stable as you can get) version will always be maintained in this repo.


### __Sections under development__

Currently the STAR and readGroup rules have not been finalised.


### Credits

Alastair Ludington
- email: alastair.ludington@adelaide.edu.au

Jimmy Breen
- email: jimmy.breen@adelaide.edu.au
