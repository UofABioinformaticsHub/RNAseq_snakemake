##---------------------------------------------------------------------------##
# RNA-Seq analysis pipeline
# @author Alastair Ludington
# @id a164524
# @email a1645424@adelaide.edu.au
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
# This is an RNA-seq analysis snakefile script for use within the Adelaide
# University Bioinformatics hub. The snakefile is adapted from Dr. Jimmy
# Breen's bash script found at the following address:
# github.com/jimmybgammyknee/gammytools/blob/master/RRI_RNAseq_pipeline.sh.
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
# Adaptations to the samplesheet column names and the tsv2yaml.sh will result
# in discrepancies between the master snakefile and config.yaml file. Be sure
# to check that the `config["..."]` commands throughout the snakefile match
# the adaptations that have been made.
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
# To run snakemake pipeline:
# Runs on python 3
# snakemake -j 2 --core 10 -s snakefile.py
    # snakemake - specifying the snakemake executable
    # -j - specifying how many jobs to run in parallel
    # --core - specifying the number of cores to use (note: -j * --core)
    # -s - snakefile
# Run specific rule: include -R followed by rule name
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##
# Development:
    # STAR:
        # Untested
    # READGROUP:
        # Not working, in development.
##---------------------------------------------------------------------------##


##---------------------------------------------------------------------------##
# Pipeline ToC:
# PART0: Targets for pipeline
# PART1: AdapterRemoval trimming/removal
# PART2: QC raw
# PART3: QC trimmed
# PART4: Hisat2/Sambamba
# PART4.5: STAR/Sambamba
# PART5: Generating readgroups
# PART6: Salmon index
# PART6: Quantification Salmon
# PART7: FeatureCounts
# PART8: Stringtie assembly/initial file formation
# PART8.5 Stringtie making merged.transcripts.gtf
# PART9: Kallisto
##---------------------------------------------------------------------------##

##---------------------------------------------------------------------------##

# Importing required utilities
import glob

## Configuration file
# Include path if not in snakefile directory
configfile: 'config.yaml'

##----------------------------------------------------------------------------##
## 0. Targets for analysis                                                    ##
##----------------------------------------------------------------------------##

## Target files
QC_raw = expand(config["BASE"] + config["QC_raw"] + '/{samples}_fastQC/', samples=config["SampleList"])
QC_trim = expand(config["BASE"] + config["QC_trim"] + '/{samples}_fastQC/', samples=config["SampleList"])
QUANT_feat = expand(config["BASE"] + config["QUANT_FC"] + '/project_genes_{samples}.txt', samples=config["SampleList"])
SAL_index = config["INDEX_DIR"]
QUANT_salmon = expand(config["BASE"] + config["QUANT_SAL"] + "/{samples}/", samples=config["SampleList"])
KAL = expand(config["BASE"] + config["KAL"] + "/{samples}_kalOut", samples=config["SampleList"])
# RG = expand(config["BAMS_RG"] + "/{samples}.bam", samples=config["SampleList"])
STRING_ass = expand(config["BASE"] + config["ASS"] + '/{samples}_assembly.gtf', samples=config["SampleList"])
STRING_geneAbund = expand(config["BASE"] + config["ASS"] + '/{samples}_gene_abund.tab', samples=config["SampleList"])
STRING_covRefs = expand(config["BASE"] + config["ASS"] + '/{samples}_cov_refs.gtf', samples=config["SampleList"])
STRING_merge = config["BASE"] + config["ASS"] + '/merged.gtf',
STRING_mergeTranscript = config["BASE"] + config["ASS"] + "/merged.transcript.gtf"

## Specify targets
rule all:
    input:
         QC_trim + QC_raw + QUANT_feat + QUANT_salmon + KAL + STRING_ass + STRING_geneAbund + STRING_covRefs, STRING_merge, STRING_mergeTranscript, SAL_index
         ## Work in progress - RG


##--------------------------------------##
## 1. Adapter removal                   ##
##--------------------------------------##

rule adpt_trimming:
    input:
        R1 = lambda wildcards: config["SampleList"][wildcards.samples]["R1"],
        R2 = lambda wildcards: config["SampleList"][wildcards.samples]["R2"]
    output:
        R1 = config["BASE"] + config["TRIM"] + "/{samples}_t1.fastq.gz",
        R2 = config["BASE"] + config["TRIM"] + "/{samples}_t2.fastq.gz"
    params:
        AdapterRemoval = config["AdapterRemoval"],
        min_qual = config["minquality_AR"],
        min_len = config["minlength_AR"],
        name = lambda wildcards: config["basename"] + "_" + config["SampleList"][wildcards.samples]["Sample"],
        excessReads = config["BASE"] + config["TRIM"] + config["DISCARD"]
    threads: config["AR_threads"]
    message:
        "AdapterRemoval - Trimming reads from fastq files"

    shell:
        """
        mkdir -p {params.excessReads}

        {params.AdapterRemoval} --file1 {input.R1} --file2 {input.R2} \
        --output1 {output.R1} --output2 {output.R2} \
        --threads {threads} --gzip --basename {params.name} --trimqualities --trimns \
        --minquality {params.min_qual} --minlength {params.min_len}

        mv *.truncated.gz {params.excessReads}
		mv *.discarded.gz {params.excessReads}
		mv *.settings {params.excessReads}
        """


##--------------------------------------##
## 2. Quality Control - raw             ##
##--------------------------------------##

rule fastqc:
    input:
        R1 = lambda wildcards: config["SampleList"][wildcards.samples]["R1"],
        R2 = lambda wildcards: config["SampleList"][wildcards.samples]["R2"]
    output:
        config["BASE"] + config["QC_raw"] + '/{samples}_fastQC/'
    params:
        file_format = config["format"],
        kmer = config["kmer"],
        out_dir = config["BASE"] + config["QC_raw"] + '/{samples}_fastQC/'
    threads: config["FQC_threads"]
    message:
        "fastQC - Quality control of raw reads"

    shell:
        """
        module load fastQC/0.11.2
        fastqc -t {threads} -k {params.kmer} --outdir {params.out_dir} -f {params.file_format} {input.R1} {input.R2}
        """


##--------------------------------------##
## 3. Quality Control - trimmed         ##
##--------------------------------------##

rule fastqc_trimmed:
    input:
        R1 = config["BASE"] + config["TRIM"] + "/{samples}_t1.fastq.gz",
        R2 = config["BASE"] + config["TRIM"] + "/{samples}_t2.fastq.gz"
    output:
        config["BASE"] + config["QC_trim"] + '/{samples}_fastQC/'
    params:
        file_format = config["format"],
        kmer = config["kmer"],
        out_dir = config["BASE"] + config["QC_trim"] + '/{samples}_fastQC/'
    threads: config["FQC_threads"]
    message:
        "fastQC - Quality control of trimmed reads"

    shell:
        """
        module load fastQC/0.11.2
        fastqc -t {threads} -k {params.kmer} -f {params.file_format} --outdir {params.out_dir}  {input.R1} {input.R2}
        """


##--------------------------------------##
## 4. Alignment                         ##
##--------------------------------------##

rule hisat2_sambamba:
    input:
        R1 = config["BASE"] + config["TRIM"] + "/{samples}_t1.fastq.gz",
        R2 = config["BASE"] + config["TRIM"] + "/{samples}_t2.fastq.gz"
    output:
        unmapped = config["BASE"] + config["BAMS"] + config["unmapped"] + '/{samples}_unmapped.fastq',
        B1 = temp(config["BASE"] + config["BAMS"] + '/{samples}_temp.bam'),
        B2 = config["BASE"] + config["BAMS"] + '/{samples}.bam'
    params:
        hisat2 = config["hisat2"],
        sambamba = config["sambamba"],
        idx = config["INDEX_hst2"],
        qual = config["quality"]
    threads: config["H2_threads"]
    message:
        "Hisat2/Sambamba - Aligning reads and generating a merged, sorted bam file"
    log:
        config["BASE"] + config["LOG"] + '/Align_hisat2/{samples}.log'

    shell:
        """
        ({params.hisat2} -p {threads} -x {params.idx} --un-gz {output.unmapped} -1 {input.R1} -2 {input.R2} | \
        samtools view -bhS -q{params.qual} - 1> {output.B1}) 2> {log}
        {params.sambamba} sort -p -t {threads} -o {output.B2} {output.B1}
        """


##--------------------------------------##
## 4.5. Alignment                         ##
##--------------------------------------##

# rule starAlignment:
#     input:
#         R1 = config["BASE"] + config["TRIM"] + "/{samples}_t1.fastq.gz",
#         R2 = config["BASE"] + config["TRIM"] + "/{samples}_t2.fastq.gz"
#     output:
#         B1 = temp(config["BASE"] + config["BAMS"] + '/{samples}_temp.bam'),
#         B2 = config["BASE"] + config["BAMS"] + '/{samples}.bam'
#     params:
#         sambamba = config["sambamba"],
#         genomeDir = config["genomeDir"],
#         compress = config["compress"]
#     threads: config["ST_threads"]
#     message: 'STAR: Alignment of trimmed reads'
#     log: config["BASE"] + config["LOG"] + '/STAR/{samples}.log'
#     shell:
#         """
#         module load STAR/2.5.1a-foss-2015b
#
#         (STAR --runThreadN {threads} \
#                 --genomeDir {params.genomeDir} \
#                 --readFilesIn {input.R1} {input.R2} \
#                 --readFilesCommand gunzip -c \
#                 --outSAMtype BAM Unsorted \
#                 --outBAMcompression {params.compress} | \
#                 samtools view -bhS -q{params.qual} - 1> {output.B1}) 2> {log}
#         {params.sambamba} sort -p -t {threads} -o {output.B2} {output.B1}
#         """


##--------------------------------------##
## 5. ReadGroup - bam headers         ##
##--------------------------------------##

rule read_group:
    input:
        config["BASE"] + config["BAMS"] + "/{samples}.bam"
    output:
        config["BASE"] + config["BAMS_RG"] + "/{samples}.bam"
    params:
        picard = config["picard"],
        RGID = lambda wildcards: config["SampleList"][
            wildcards.samples]["Replicate"],
        RGLB = lambda wildcards: config["SampleList"][
            wildcards.samples]["Treatment"],
        RGPL = lambda wildcards: config["SampleList"][
            wildcards.samples]["Platform"],
        RGPU = lambda wildcards: config[
            "SampleList"][wildcards.samples]["Index"],
        RGSM = lambda wildcards: config[
            "SampleList"][wildcards.samples]["Sample"]
    # threads: config["threads"]
    message:
        "Picard AddOrReplaceReadGroups - Generating bam file readgroups"

    shell:
        """
        java -jar {params.picard} AddOrReplaceReadGroups I={input} O={output} RGID={params.RGID} RGLB={params.RGLB} RGPL={params.RGPU} RGPU={params.RGPU} RGSM={params.RGSM}
        """

##--------------------------------------##
## 6. Quantification - Salmon Index     ##
##--------------------------------------##

rule salmon_index:
    input:
        Tran = config["TRAN"]
    output:
        index_out = config["INDEX_DIR"]
    params:
        salmon = config["salmon"],
        index_out = config["INDEX_DIR"]
    threads: 1
    message: "Salmon - Generating index from transcripts file"
    shell:
        """
        {params.salmon} index -p {threads} -i {params.index_out} -t {input.Tran} --gencode --perfectHash
        """

##--------------------------------------##
## 6. Quantification - Salmon           ##
##--------------------------------------##

rule quant_salmon:
    input:
        R1 = config["BASE"] + config["TRIM"] + "/{samples}_t1.fastq.gz",
        R2 = config["BASE"] + config["TRIM"] + "/{samples}_t2.fastq.gz",
        index_in = rules.salmon_index.output
    output:
        out = config["BASE"] + config["QUANT_SAL"] + "/{samples}/"
    params:
        salmon = config["salmon"],
        bootstrap = config["bootstrap"],
        out = config["BASE"] + config["QUANT_SAL"] + "/{samples}/"
    threads: config["SAL_threads"]
    message:
        "Salmon - Quantification of transcripts"
    log:
        config["BASE"] + config["LOG"] + config["QUANT_SAL"] + "/{samples}.log"

    shell:
        """
        ({params.salmon} quant -l A -i {input.index_in} -p {threads} -1 {input.R1} -2 {input.R2} --numBootstraps {params.bootstrap} -o {params.out}) 2> {log}
        """


##--------------------------------------##
## 7. Quantification - FeatureCounts    ##
##--------------------------------------##

rule quant_featCount:
    input:
        B1 = config["BASE"] + config["BAMS"] + "/{samples}.bam",
        gtf = config["GTF"]
    output:
        TPG = config["BASE"] + config["QUANT_FC"] + '/project_genes_{samples}.out',
        PG = config["BASE"] + config["QUANT_FC"] + '/project_genes_{samples}.txt'
    params:
        featurecounts = config["featurecounts"],
        qual = config["minquality_FC"],
        strandness = config["strandness"]
    threads: config["FC_threads"]
    message:
        "Featurecounts - Transcript quantification at the gene level"
    log:
        config["BASE"] + config["LOG"] + '/Quant_featureCount/{samples}.log'

    shell:
        """
        ({params.featurecounts} -Q {params.qual} -s {params.strandness} -T {threads} -a {input.gtf} -o {output.TPG} {input.B1}
        cut -f1,7- {output.TPG} | sed 1d > {output.PG}) 2> {log}
        """


##--------------------------------------##
## 8. Assembly - Part 1                 ##
##--------------------------------------##

rule assembly_stringtie:
    input:
        BAM = config["BASE"] + config["BAMS"] + '/{samples}.bam',
        gtf = config["GTF"]
    output:
        ass = config["BASE"] + config["ASS"] + '/{samples}_assembly.gtf',
        genabund = config["BASE"] + config["ASS"] + '/{samples}_gene_abund.tab',
        covref = config["BASE"] + config["ASS"] + '/{samples}_cov_refs.gtf'
    params:
        stringtie = config["stringtie"],
        tran_len = config["minlength_ST"]
    threads: config["ST_threads"]
    message:
        "Stringtie pt 1 - Generating assembly, gene abundance and coverage statistics"

    shell:
        """
        {params.stringtie} -p {threads} -m {params.tran_len} -G {input.gtf} \
        -o {output.ass} \
        -A {output.genabund} \
        -C {output.covref} \
        {input.BAM}
        """


##--------------------------------------##
## 8.5 Stringtie merging assembly files ##
##--------------------------------------##

## Gets all files that have the *_assembly.gtf string
files = glob.glob(config["BASE"] + config["ASS"] + "/*_assembly.gtf")

rule stringtie_merge:
    input:
        gtfFiles = files,
        gtf = config["GTF"]
    output:
        mergeGTF = config["BASE"] + config["ASS"] + '/merged.gtf',
        mergedTRAN = config["BASE"] + config["ASS"] + "/merged.transcript.gtf"
    params:
        stringtie = config["stringtie"]
    message:
        "Stringtie pt 2 - Generating a merged gtf file from all samples"
    shell:
        """
        {params.stringtie} --merge -G {input.gtf} -o {output.mergeGTF} {input.gtfFiles}
        awk '{{if($3=="transcript")print}}' {output.mergeGTF} > {output.mergedTRAN}
        """


##--------------------------------------##
## 9. Kallisto                          ##
##--------------------------------------##

rule kallisto:
    input:
        R1 = config["BASE"] + config["TRIM"] + "/{samples}_t1.fastq.gz",
        R2 = config["BASE"] + config["TRIM"] + "/{samples}_t2.fastq.gz",
        KAL = config["KAL_idx"]
    output:
        out = config["BASE"] + config["KAL"] + "/{samples}_kalOut"
    params:
        kallisto = config["kallisto"],
        bootstrap = config["bootstrap"],
        outDir = config["BASE"] + config["KAL"] + "/{samples}_kalOut"
    threads: config["KAL_threads"]
    message: 'Kallisto - transcript quantification (validation)'

    shell:
        """
        {params.kallisto} quant -i {input.KAL} -t {threads} -b {params.bootstrap} -o {params.outDir} {input.R1} {input.R2}
        """
