###########################################################################
## Author: Alastair Ludington
## Email: alastair.ludington@adelaide.edu.au
## Date: 15/01/2018
## Developed on behalf of the Bioinformatics Hub, University of Adelaide
##
## The purpose of this pipeline is to generate gene count data from a
## range of quantification software.
## The pipeline is still in active development. This is v0.1, which is 
## usable but still has growing pains
###########################################################################

# from snakemake.utils import R
# import glob

configfile: 'config.yaml'

QC_raw = config["BASE"] + config["FASTQC"] + config["QC_raw"]
QC_trim = config["BASE"] + config["FASTQC"] + config["QC_trim"]
SAL = expand(config["BASE"] + config["QUANT_SAL"] + '/{samples}/', samples=config["SampleList"])
KAL = expand(config["BASE"] + config["QUANT_KAL"] + "/{samples}/", samples=config["SampleList"])
FC = config["BASE"] + config["QUANT_FC"] + '/counts.txt'
ST = config["BASE"] + config["QUANT_STR"] + "/merged.transcript.gtf"
STAR = expand(config["BASE"] + config["ALIGN"] + config["BAMS"] + '/{samples}.bam', samples=config["SampleList"])
# HTML = config["BASE"] + config["FASTQC"] + config["QC_trim"] + '/ngsReports_Fastqc.html'


## add/remove variables that you want to run
rule all:
    input:
    	STAR + SAL + KAL, FC, ST, QC_raw, QC_trim
        


##--------------------------------------##
## 1. Adapter removal                   ##
##--------------------------------------##

rule Adapter_removal:
    input:
        R1 = lambda wildcards: config["SampleList"][wildcards.samples]["R1"],
        R2 = lambda wildcards: config["SampleList"][wildcards.samples]["R2"]
    output:
        R1 = config["BASE"] + config["AR"] + "/{samples}_trim1.fastq.gz",
        R2 = config["BASE"] + config["AR"] + "/{samples}_trim2.fastq.gz"
    params:
        min_qual = config["minquality_AR"],
        min_len = config["minlength_AR"],
    threads: config["THREADS"]
    message:
        "ADAPTERREMOVAL - Trimming reads from file {wildcards.samples}"

    shell:
        """
        AdapterRemoval \
        --file1 {input.R1} \
        --file2 {input.R2} \
        --output1 {output.R1} \
        --output2 {output.R2} \
        --threads {threads} \
        --gzip \
        --trimqualities \
        --trimns \
        --minquality {params.min_qual} \
        --minlength {params.min_len}
        """


##--------------------------------------##
## 2. Quality Control - raw             ##
##--------------------------------------##

SUF = ["_R1","_R2"]

rule fastqc_raw:
    input:
        expand(config["RAW_FILE_PATH"] + '/{samples}{suf}.fastq.gz', samples=config["SampleList"], suf=SUF)
    output:
        config["BASE"] + config["FASTQC"] + config["QC_raw"]
    params:
        file_format = config["format"],
        kmer = config["kmer"]
        # RMD = config["BASE"] + '/ngsReports_Fastqc.Rmd'
    threads: config["THREADS"]
    message:
        "FASTQC - Quality control of raw reads"

    shell:
        """
        mkdir -p {output}
        fastqc \
        -t {threads} \
        -k {params.kmer} \
        --outdir {output} \
        -f {params.file_format} \
        {input}
        """
#       cp {params.RMD} {output}

##--------------------------------------##
## 3. Quality Control - trimmed         ##
##--------------------------------------##

TRIM  = ["trim1", "trim2"]

rule fastqc_trimmed:
    input:
        expand(config["BASE"] + config["AR"] + '/{samples}_{trim}.fastq.gz', samples=config["SampleList"], trim=TRIM)
    output:
        config["BASE"] + config["FASTQC"] + config["QC_trim"]
    params:
        file_format = config["format"],
        kmer = config["kmer"]
        # RMD = config["BASE"] + '/ngsReports_Fastqc.Rmd'
    threads: config["THREADS"]

    shell:
        """
        mkdir -p {output}

        fastqc \
        -t {threads} \
        -k {params.kmer} \
        -f {params.file_format} \
        --outdir {output} \
        {input}
        """
#        cp {params.RMD} {output}

# rule FQC_report:
#     input: 
#         RMD = config["BASE"] + config["FASTQC"] + config["QC_trim"] + '/ngsReports_Fastqc.Rmd'
#     output:
#         config["BASE"] + config["FASTQC"] + config["QC_trim"] + 'ngsReports_Fastqc.html'
#     params:
#         NAME = 'ngsReports_Fastqc.html',
#         PATH = config["BASE"] + config["FASTQC"] + config["QC_trim"]
#     run:
#         R("""
#         library(rmarkdown)
#         render({input})  
#         """)

# rmarkdown::render({input}, output_format = "html_document",
#         output_file = {params.NAME},output_dir = {params.PATH},
#         knit_root_dir = dirname({input}), envir = new.env(), quiet = quiet,
#         params = list(dataType = "Transcriptome", species = "Hsapiens"))


##--------------------------------------##
## 4. STAR ALIGNMENT                    ##
##--------------------------------------##

rule STAR_index:
    input:
        genomeFasta = config["TRAN"],
        GTFfile = config["GTF"]
    output:
        config["BASE"] + config["ALIGN"] + '/STAR_index'
    params:
        overhang = config["overhang"]
    threads: config["THREADS"]
    message: 
        "STAR - Generating Index for {input.genomeFasta} using {input.GTFfile}"
    log:
        config["BASE"] + config["LOG"] + '/STAR_index.log'
    shell:
        """
        STAR \
        --runMode genomeGenerate \
        --runThreadN {threads} \
        --genomeDir {output} \
        --genomeFastaFiles {input.genomeFasta} \
        --sjdbGTFfile {input.GTFfile} \
        --sjdbOverhang {params.overhang} 2> {log}
        """


rule STAR_alignment:
    input:
        R1 = config["BASE"] + config["AR"] + "/{samples}_trim1.fastq.gz",
        R2 = config["BASE"] + config["AR"] + "/{samples}_trim2.fastq.gz",
        IDX = rules.STAR_index.output
    output:
        B0 = config["BASE"] + config["ALIGN"] + config["BAMS"] + '/{samples}_STAR.Aligned.out.bam',
        B1 = config["BASE"] + config["ALIGN"] + config["BAMS"] + '/{samples}_filtered.bam',
        B2 = config["BASE"] + config["ALIGN"] + config["BAMS"] + '/{samples}.bam'
    params:
        BAMcompression = config["BAMcompression"],
        qualityFilter= config["qualityFilter"],
        OUT = config["BASE"] + config["ALIGN"] + config["BAMS"] + '/{samples}_STAR.'
    threads: config["THREADS"]
    message: 
    	'STAR - Aligning reads of {wildcards.samples}'
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --genomeDir {input.IDX} \
        --readFilesIn {input.R1} {input.R2} \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix {params.OUT} \
        --outSAMtype BAM Unsorted \
        --outBAMcompression {params.BAMcompression} 

        samtools view -h -b -q {params.qualityFilter} {output.B0} -o {output.B1}
        
        sambamba sort -p -t {threads} -o {output.B2} {output.B1}

        """

        
##--------------------------------------##
## 5. ReadGroup - bam headers         	##
##--------------------------------------##

rule read_group:
    input:
        config["BASE"] + config["ALIGN"] + config["BAMS"] + '/{samples}.bam'
    output:
        RG = config["BASE"] + config["ALIGN"] + config["RG_PICARD"] + '/{samples}.bam',
        IDX = config["BASE"] + config["ALIGN"] + config["RG_PICARD"] + '/{samples}.bam.bai'
    params:
        RGID = lambda wildcards: config["SampleList"][wildcards.samples]["RG_ID"],
        RGLB = lambda wildcards: config["SampleList"][wildcards.samples]["RG_LB"],
        RGPL = lambda wildcards: config["SampleList"][wildcards.samples]["RG_PL"],
        RGPU = lambda wildcards: config["SampleList"][wildcards.samples]["RG_PU"],
        RGSM = lambda wildcards: config["SampleList"][wildcards.samples]["RG_SM"]
    message:
        "PICARD - Assigning readgroup information to {wildcards.samples}"
    threads: config["THREADS"]
    log:
        config["BASE"] + config["LOG"] + config["RG_PICARD"] + '/{samples}.log'

    shell:
        """
        java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
        I={input} \
        O={output.RG} \
        RGID={params.RGID} \
        RGLB={params.RGLB} \
        RGPL={params.RGPL} \
        RGPU={params.RGPU} \
        RGSM={params.RGSM}

        sambamba index -t {threads} {output.RG} {output.IDX}
        """


##--------------------------------------##
## 6. Salmon - Quantification        	##
##--------------------------------------##

rule salmon_index:
	input:
		config["TRAN"]
	output:
		config["BASE"] + config["QUANT_SAL"] + config["SALMON_INDEX"]
	threads: config["THREADS"]
	message: 
		"SALMON - Generating Index of {input} and storing it in {output}"
	shell:
	    """
	    salmon index \
	    -p {threads} \
	    -i {output} \
	    -t {input} \
	    --gencode \
	    --perfectHash 
	    """

rule quant_salmon:
    input:
        R1 = config["BASE"] + config["AR"] + "/{samples}_trim1.fastq.gz",
        R2 = config["BASE"] + config["AR"] + "/{samples}_trim2.fastq.gz",
        index = rules.salmon_index.output
    output:
        out = config["BASE"] + config["QUANT_SAL"] + '/{samples}/'
    params:
        bootstrap = config["bootstrap"]
    threads: config["THREADS"]
    message:
        "SALMON - Quantification of transcripts for samples {wildcards.samples}"
    log:
        config["BASE"] + config["LOG"] + config["QUANT_SAL"] + "/{samples}.log"

    shell:
        """
        salmon quant \
        -l A \
        -i {input.index} \
        -p {threads} \
        -1 {input.R1} \
        -2 {input.R2} \
        --numBootstraps {params.bootstrap} \
        -o {output} 2> {log}
        """

##--------------------------------------##
## 7. Kallisto - quantification         ##
##--------------------------------------##

rule kallisto_index:
	input:
		fasta = config["TRAN"],
	output:
		config["BASE"] + config["QUANT_KAL"] + config["KALLISTO_INDEX"]
	message:
		"KALLISTO - Generating Kallisto index file"
	shell:
		"""
        module load kallisto/0.43.1-foss-2017a
		kallisto index -i {output} \
		{input.fasta}
		"""

rule quant_kallisto:
    input:
        R1 = config["BASE"] + config["AR"] + "/{samples}_trim1.fastq.gz",
        R2 = config["BASE"] + config["AR"] + "/{samples}_trim2.fastq.gz",
        IDX = rules.kallisto_index.output
    output:
        out = config["BASE"] + config["QUANT_KAL"] + "/{samples}/"
    params:
        bootstrap = config["bootstrap"],
    threads: config["THREADS"]
    message: 'Kallisto - transcript quantification (validation)'
    log:
        config["BASE"] + config["LOG"] + config["QUANT_KAL"] + '/{samples}.log'

    shell:
        """
        module load kallisto/0.43.1-foss-2017a
        kallisto quant \
        -i {input.IDX} \
        -t {threads} \
        -b {params.bootstrap} \
        -o {output} {input.R1} {input.R2}
        """

##--------------------------------------##
## 8. FeatureCounts - Quantification    ##
##--------------------------------------##

rule quant_featCount:
    input:
    	BAMS = expand(config["BASE"] + config["ALIGN"] + config["RG_PICARD"] + '/{samples}.bam', samples=config["SampleList"]),
        GTFfile = config["GTF"]
    output:
        TPG = config["BASE"] + config["QUANT_FC"] + '/counts.out',
        PG = config["BASE"] + config["QUANT_FC"] + '/counts.txt'
    params:
        qual = config["minquality_FC"],
        strandness = config["strandness"]
    threads: config["THREADS"]
    message:
        "FEATURECOUNTS - Transcript quantification at the gene level"

    shell:
        """
        featureCounts \
        -Q {params.qual} \
        -s {params.strandness} \
        -T {threads} \
        -a {input.GTFfile} \
        -o {output.TPG} {input.BAMS}

        cut -f1,7- {output.TPG} | \
        sed 1d > {output.PG}
        """

##--------------------------------------##
## 9. Stringtie - Quantification        ##
##--------------------------------------##

rule quant_stringtie_assembly:
    input:
        BAM = config["BASE"] + config["ALIGN"] + config["RG_PICARD"] + '/{samples}.bam',
        GTFfile = config["GTF"]
    output:
        STR_OUT = config["BASE"] + config["QUANT_STR"] + '/{samples}_assembly.gtf',
        genabund = config["BASE"] + config["QUANT_STR"] + '/{samples}.tab',
        covref = config["BASE"] + config["QUANT_STR"] + '/{samples}_cov_refs.gtf'
    params:
        tran_len = config["minlength_ST"]
    threads: config["THREADS"]
    message:
        "STRINGTIE  - Gene abundance and coverage statistics for {wildcards.samples}"

    shell:
        """
        stringtie -p {threads} \
        -m {params.tran_len} \
        -G {input.GTFfile} \
        -o {output.STR_OUT} \
        -A {output.genabund} \
        -C {output.covref} \
        {input.BAM}
        """


rule quant_stringtie_merge:
    input:
        ST_GTF = expand(config["BASE"] + config["QUANT_STR"] + '/{samples}_assembly.gtf', samples=config["SampleList"]),
        GFTfile = config["GTF"]
    output:
        mergeGTF = config["BASE"] + config["QUANT_STR"] + '/merged.gtf',
        mergedTRAN = config["BASE"] + config["QUANT_STR"] + "/merged.transcript.gtf"
    message:
        "STRINGTIE - Generating merged GTF from all samples"
    shell:
        """
        stringtie --merge -G {input.GFTfile} -o {output.mergeGTF} {input.ST_GTF}
        awk '{{if($3=="transcript")print}}' {output.mergeGTF} > {output.mergedTRAN}
        """