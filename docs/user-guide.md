# Project Name User Guide v0.0.0

## Introduction
Put an introduction here.

## Main Content
Put more content here.


  

# iCOMIC

iCOMIC (**I**ntegrating **C**ontext **O**f **M**utation **I**n **C**ancer) an open-source, standalone tool for genomic data analysis characterized by a Python-based Graphical User Interface and automated Bioinformatics pipelines for DNA-Seq and RNA-Seq data analysis. It serves as a point and click application facilitating genomic data analysis accessible to researchers with minimal programming expertise.

## Content

  

 1. About 
	1. Major features
 2. Installation
	1. Prerequisites
	2. Github download and installation
	3. Requirements.txt
	4. Conda Installation 
 3. Getting Started
	1. iCOMIC overview
	2. Install iCOMIC
	3. Usage
	4. Terminology
	5. Launching the wrapper
	6. Input file format
	7. Quick Guide
	8. Output Information
	9. FAQs
 4. Using iCOMIC for DNA-Seq or RNA-Seq analysis
	1. Analysis steps in DNA Seq
		1. Input Requirements
		2. Review of input samples
		3. Setting up your custom pipeline
		4. Initialization of the analysis
		5. Results a quick check
	2. Analysis steps in RNA Seq 
		1. Input Requirements
		2. Review of input samples
		3. Setting up your custom pipeline
		4. Initialization of the analysis
		5. Results a quick check

 5. Walkthrough of pre-constructed pipeline
	1. List of pipelines
	2. Description of the tools used 

 6. Creating a custom pipeline

	 1. Shell scripts to be written by the user
	 2. How to run a custom pipeline

 7. Viewing and analyzing results
	 1.  MultiQC Reports
	 2. Plots generated in RNA Seq
	 3. List of differentially expressed genes
	 4. Variants in DNA Seq
	 5. Annotated variants in DNA Seq

 8. Gallery
	 1. Screenshots
	 2. GIF

 9. Case examples
	 1. Installation
	 2. Configuring iCOMIC
	 3. Analysis using test data WGS
	 4. Analysis using test data RNA Seq
	 5. Results 

 10. Troubleshooting runtime issues
	 1. FAQs

 11. Changelog

 12. Glossary
## About
iCOMIC is a user-friendly pipeline to analyze genomic data that takes in raw sequencing data in FASTQ format as an input, and finally outputs insightful statistics on it’s nature. iCOMIC toolkit is capable of analyzing both Whole Genome and transcriptome data and is embedded in ‘Snakemake’, a workflow management system. iCOMIC is characterized by a user-friendly GUI built using PyQt5 which increases its ease of access. The toolkit features many independent core workflows in both whole genomic and transcriptomic data analysis pipelines.

#### Major features

 -    Serves as a stand-alone end to end analysis toolkit for DNA-Seq and RNA-Seq data
    
-   Characterized by an interactive and user friendly GUI, specifically built to accommodate users with minimal programming expertise
    
-   Provides expert bioinformaticians a platform to perform analysis incorporating advanced parameters, saving time on building a pipeline
    
-   Consists of multiple flexible workflows
    
-   Users freedom to select tools from the predesigned combinations best suited for their requirements
    
-   Easy installation of tools and dependencies

## Installation

  

Setting up iCOMIC is comparatively effortless across Linux or Mac platforms. iCOMIC works with Python 3.6 or above.

#### 1. Prerequisites

- Linux/Windows/Mac platform

- Python 3.6 and above

- iCOMIC package downloaded from GitHub

- Memory requirement

#### 2. Github download and installation

The entire source code for the tool is available at “ ”

#### 3. Requirements.txt

Install the dependencies associated with iCOMIC by using the following command

```

$ pip install requirements.txt

```

  #### 4. Conda Installation

## Getting started

This guide will walk you through the steps necessary to understand, install, and use iCOMIC for carrying out analysis on your data.

  

#### 1. iCOMIC overview

iCOMIC is an open-source, stand-alone toolkit for genomic data analysis, characterized by a python based Graphical User Interface. The tool enables researchers with minimal programming expertise to draw consequential insights from DNA Seq and RNA Seq data.

#### 2. Install iCOMIC

Installation is easy as we provide a `requirements.txt` file comprising all the software dependencies. Once you clone the iCOMIC github repository, you can install all the associated dependencies using the command below. Every additional software requirement will be managed by the conda environment.

```

$ pip install requirements.txt

```

  

#### 4. Usage

1) Snakemake

	iCOMIC is embedded in Snakemake, a Python based workflow manager. Different tools integrated in iCOMIC are connected using Snakemake. Individual ‘Rules’ corresponding to each tool form the building units, which describes how the desired output is obtained from the input. Rules consist of information about the input and output files and wrapper script or shell command. Tools without wrapper scripts are configured separately and shell command is used for their execution. According to the choice of tools made by the user, corresponding rules are combined in a ‘Snakefile’ to generate target output. All the input information and parameters corresponding to each tool is specified in a configuration file, ‘config file’. ‘Rules’ are predefined and are made available together with the iCOMIC package. All the other files are generated on the flow according to the user inputs and are updated accordingly. 

2) PyQt5 GUI

	iCOMIC is characterized by a Graphical user Interface which enhances the accessibility of the toolkit. The GUI framework is built using PyQt5, python binding of the cross-platform GUI toolkit Qt. The GUI framework allows users with minimal programming expertise to perform analysis.

  

#### 5. Terminology

  

#### 6. Launching the wrapper

iCOMIC can be launched using a simple command in the terminal.

```

$ ‘python mainwin_v32.py’

```

#### 7. Input file format

iCOMIC accepts input information in two different modes. The user can either feed the path to a folder containing raw fastq files or provide a table consolidating particulars of raw data.

  

>If you are uploading a folder of fastq files, all the files should be named in the specified format:

(sample_name)_(condition(tumor/normal))_Rep(replicate_number)_R(1 / 2).fastq

Example: hcc1395_normal_Rep1_R1.fastq

  

>If you choose upload from table mode, the sample information should be given in a tab delimited file with a header row.

The Column names should be:

`Sample` : The Sample name

`Unit` : The number of replicates

`Condition` : Nature of the sample, Normal or Tumor

`fq1` : The path of Read 1

`fq2` : Path of Read 2, if you are working with single-end reads only, the 'fq2' column can be left blank.

  

#### 8. Quick Guide

iCOMIC toolkit enables the ready analysis of RNA-Seq and Whole Genome Sequencing data. It has an inbuilt library of tools with predefined valid combinations. The user will have the freedom to choose any possible combination of tools. Figure 1 and 2 depicts the basic steps and outputs involved in DNA Seq and RNA Seq pipelines respectively.
![Figure 1: DNA Seq pipeline](https://lh3.googleusercontent.com/N1xl3xEPuLhDBvQucuzAEoAoVnzktLMyZdiDloSNIo0UIxE3k0B_aAwf3PcKdmZNDZ4O8zCCQA3_W22VCGc7P63xXWFRbYhE5ESRq9t5dQ63xR_Glh5nRhPHaDJiHJHIwtJP26w1p1gJ9wf8lXy3OkaYakuD5XHlfsqLo306Cp8O14UL8Ng4oKYuaw2XJyGK4q9LNg6waSXlcirD_3G0QUsz4ufhqcahS8opW-p2NdA7z_y7moB_-UDs81pjCbpq2rL2sm5R211WwNPtue5Nbz353j6LR3ZWjweGRx6wui5tRNE7_ei9SSKEAtE7fOM6oa6zh-b_XfXcPdKFu805UFCvvbaO321xL2G-TelfqDUkOcn2NbX41_LVzlAoelbwcvNvPgQ5deWoO2c_VOgIBa58cQzruM1-fNDeC7uuaHoTyAKV2CDTuH8oU3KBKSBu1MmwJFuojoseBsgFNhq52dO7CqzG8EgMvi4T0ANzw0kx4z7LBQZyLXAOy-PYexvVfHP3789e_aSjkMyZQp8XCbQr7-6Js4Tfq9D6D763pn23s4Ekf63aplEt0g8BY-LbnmkJgEN704qMTCtg-t3GB6A9BkM7rTQXpshBF1i4tPXc4eYFcnRR6NPEuwZvsLYsrY2hsj_Qncq5dztPHbmq-p6fwP_0yMzGpgrzdu7GWAzQw3l68OdSUOkOVzRfm6I=w1666-h937-no)  
Figure 1: DNA Seq pipeline. This workflow indicates the analysis steps and major output files in DNA Seq pipeline.
![enter image description here](https://lh3.googleusercontent.com/l-Lp8JkEP3-MkFqsdQkj6TIYW0kuVTDrngCLLgtzgQnnEeb9XLLlwyArCnyDx8tsawLRm8q8kiBbrLtt7EegchKCX4jQ_8rJfIhr5nfj9Xhq6naqUZIZxDelTRev3yadGV5ujChPVIwYkFM4FWq1ICcbln072yNVGSGjkZGNoFTwDICnmpLwkwwEEiQh4OwlfLDMqepamVM4GueUU3VYAn9Lhea7FjuqB90HqDpLdsHARgy9lZ1uWBMPWG4HS1Bz_SSmvk09XM8uToTs7n9FCp6DpZvd63qmZPHteHmTjfVP8BfmwTtf17JucG8rZSUoENjvZhx-0faNvtDFP3KtiDnAHyV57h1_cP2LNGZQSzHX9iWCCuyGXLZ1yv0veYn5P3ai-10EwwKVQyC3wae2R82h_pyCBOloy99v2Gwla1rXOa438lbsftEVLQZDNEY0i8t-gkVrFoEnAnT00QYGM6PHGcdUy0JgR7WD1YwjGg8FHXlreDeiYxaUP9L5TyjhHVQyfU4NOO4oP3XHsHzxaFHwOi5AkJuttmTEG2U1IrOmdyDVlNkOElHmA0Kw9H59jJUe-G_cNTmrEESFhsgnBLAJzIgw_ky5iElh1A3dCX2CNw4KO_s9GVAYBuqYmNPuH8rG6r4NLbU0Wm59neKRMhv7ny1xlZmY6G34w1w85T20JmblmGZs0cDSLaWSO7U=w1666-h937-no)
Figure 2: The analysis steps and major output files in RNA Seq pipeline.
Here is a typical set of actions to run iCOMIC pipelines:
- Select a pipeline.

- Choose the mode of input

- Input the required data fields.

- Proceed to the next tab if you want to skip Quality Check.

- Or click on the `Quality Control Results` button to view a consolidated MultiQC report of Quality statistics.

- Check `yes` if you want to do trimming and also mention the additional parameters as per requirement.

- Tool for Quality Control: FastQC

- Tool for trimming the reads: Cutadapt

- Choose the tools of interest from `Tool selection` tab and set the parameters as required

- For the choice of aligner, the corresponding genome index file needs to be uploaded if available, or the user can generate the index file using the `Generate Index` button.

- Click `Run` on the next tab to run the analysis.

- If a warning button pops-up near the `Unlock` button, click on it to unlock the working directory.

- Once the analysis is completed, `Results` tab will be opened.

- DNA Seq results include a MultiQC report comprising the statistics of the entire analysis. A file consisting of the variants called and the corresponding annotated variant file.

- Results for RNA Seq analysis include multiQC analysis statistics, R plots such as MA plot, Heatmap, PCA plot and box plot and list of differentially expressed genes.

#### 9. Output information

All outputs are stored in separate folders for each pipeline along with log information.

- DNA-Seq

	- MultiQC
	Contains subfolders MultiQC, FastQC and Cutadapt. MultiQC contains consolidated `html` reports on the overall run statistics and a separate `html` file on merged FastQC reports of all the input samples. The folder FastQC contains quality reports of individual samples. It may also enclose FastQC report of trimmed reads if the user opts for trimming the input reads. The folder Cutadapt contains trimmed `fastq` files.
	- Aligner
Contain `bam` outputs generated by the aligner.
	- Variant Caller
This folder includes `vcf` files of identified variants.
	- Annotator
Contains annotated `vcf` files.
	- Index
This is an optional folder which contains the index files if the user chooses to generate index corresponding to the choice of aligner.
- RNA-Seq

	- MultiQC
Contains subfolders MultiQC, FastQC and Cutadapt. MultiQC contains consolidated `html` reports on the overall run statistics and a separate `html` file on merged FastQC reports of all the input samples. The folder FastQC contains quality reports of individual samples. It may also enclose FastQC report of trimmed reads if the user opts for trimming the input reads. The folder Cutadapt contains trimmed `fastq` files.
	- Aligner
Contain `bam` outputs generated by the aligner.
	- Expression Modeller
Contains count matrix representing the reads mapped to individual genes.
	- Differential expression
Contain a text file with a consolidated list of differentially expressed genes.
	- Index
This is an optional folder which contains the index files if the user chooses to generate index corresponding to the choice of aligner.
#### 10. FAQs

- What to do if run fails

- How to run iCOMIC?

- Dependences

- Installation issues

  

## Using iCOMIC for DNA-Seq or RNA-Seq analysis

#### 1. Analysis steps DNA Seq

`DNA Seq` constitutes the Whole Genome Sequencing pipeline which permits the user to call variants from the input samples and annotate them. iCOMIC integrates a combination of 3 aligners, 5 variant callers and 2 annotators along with the tools for Quality control. The tool MultiQC is incorporated to render comprehensive analysis statistics.

  

[Can give a sample DAG]

  

##### 1. Input Requirements

The significant obligation is raw `fastq` files which can either be single-end or paired-end. Fastq read details can be specified in two different methods, either by uploading a folder containing the reads or using a tab-separated file describing the reads. If you choose the `Upload from folder` mode, the path to the folder containing the fastq files needs to be specified. Additionally, all the files in the folder need to be named in the format specified in the section ` ‘Input file format’ `.

Alternatively, if you decide to use the ‘Upload from table’ mode, a tab-separated file consolidating particulars about the sample needs to be fed in. Refer to the ` ‘Input file format’` section for formatting the table.

  

Furthermore, iCOMIC demands a path to the reference genome, a `fasta` file and gzipped `vcf` of the known variants corresponding to the reference genome.

  

Once all the fields are filled, you can proceed to the `Quality Control` tab using the next button.

##### 2. Review of Input Samples

In the `Quality Control` widget, you can examine the quality of your samples for analysis by clicking on the `Quality Control Results` button. The tool MultiQC provides a consolidated report of Quality statistics generated by FastQC for all the samples. Additionally, iCOMIC permits you to trim the reads using Cutadapt if required. However one will be free to move ahead even without going through the Quality Check processes.

##### 3. Setting up your custom pipeline

In the tool selection widget, you will be asked to choose your desired set of tools for analysis.

- Aligner

You can choose a software for sequence alignment from the drop down menu. You will also need to input the genome index corresponding to the choice of aligner. No worries! iCOMIC allows you to generate the required index using the `Generate index` button. One will have the permission to change the values for the mandatory parameters displayed. Moreover, if you are an expert bioinformatician, iCOMIC allows you to play around with the advanced parameters. Clicking on the `Advanced` button would open a pop-up of all the parameters associated with a tool.

- Variant Caller

This section permits you to choose a variant caller from the set of tools integrated. If the input sample is normal-tumor specific, then only those tools which call variants comparing the normal and tumor samples will be displayed. On the other hand, if you want to call variants corresponding to the reference genome, variant callers of that type would be displayed. iCOMIC allows you to set mandatory as well as advanced parameters for the selected tool.

- Annotator

This section allows you to choose a tool for annotating your called variants and specify the parameters.

  

##### 4. Initialization of the analysis

The `Run` tab displays an `Unlock` button and a `Run` button. `Run` is for initializing the analysis. When the analysis starts, if a warning icon pops up near the `Unlock` button, you need to click the `Unlock` button to unlock the working directory and then click `Run` to proceed with the analysis. Progress bar present in the tab allows you to examine the progress of analysis.

##### 5. Results: a quick check

Once the analysis is completed, iCOMIC will automatically take you to the `Results` tab which displays three major results.

1. Analysis Statistics

Displays a MultiQC consolidated report of overall analysis statistics. This includes FastQC reports, Alignment statistics and variant statistics.

2. Variants called

On clicking this button a pop up with the `vcf` file of variants called will be displayed.

3. Annotated variants

Displays the annotated `vcf` file.

  

#### 2. Analysis steps RNA Seq

`RNA Seq` part allows you to identify the differentially expressed genes from RNA Sequencing data.iCOMIC integrates a combination of 3 aligners, 3 expression modellers and 2 differential expression tools along with the tools for Quality control. The tool MultiQC is incorporated to render comprehensive analysis statistics.

  

[Can give a sample DAG]

  

##### 1. Input Requirements

Similar to the DNA Seq pipeline, the major requirement here is also raw `fastq` files, either single-end or paired-end. Fastq read details can be specified in two different methods, either by uploading a folder containing the reads or using a tab-separated file describing the reads. Refer to the `‘Input file format’` section for preparing the input reads. Furthermore, the RNA Seq part demands a path to the reference genome, a `fasta` file, annotated file in `gtf` format, and a transcript file.

Once all the fields are filled, you can proceed to the `Quality Control` tab using the next button.

  

##### 2. Review of Input Samples

In the `Quality Control` widget, you can examine the quality of your samples for analysis by clicking on the `Quality Control Results` button. The tool MultiQC provides a consolidated report of Quality statistics generated by FastQC for all the samples. Additionally, iCOMIC permits you to trim the reads using Cutadapt if required. However one will be free to move ahead even without going through the Quality Check processes.

  

##### 3. Setting up your custom pipeline

In the tool selection widget, you will be asked to choose your desired set of tools for analysis.

- Aligner

You can choose a software for sequence alignment from the drop down menu. You will also need to input the genome index corresponding to the choice of aligner. No worries! iCOMIC allows you to generate the required index using the `Generate index` button. One will have the permission to change the values for the mandatory parameters displayed. Moreover, if you are an expert bioinformatician, iCOMIC allows you to play around with the advanced parameters. Clicking on the `Advanced` button would open a pop-up of all the parameters associated with a tool.

- Expression Modeller

This section allows you to choose an expression modeller from the integrated list of tools for counting the reads with the help of annotation file. Users will also have the freedom to set parameters corresponding to the tool.

- Differential Expression tool

Here you can choose a tool for quantifying differential expression and can also set parameters.

##### 4. Initialization of the analysis

The `Run` tab displays an `Unlock` button and a `Run` button. `Run` is for initializing the analysis. When the analysis starts, if a warning icon pops up near the unlock button, you need to click the unlock button to unlock the working directory and then click run to proceed with the analysis. Progress bar present in the tab allows you to examine the progress of analysis.

  

##### 5. Results: a quick check

Once the analysis is completed, iCOMIC will automatically take you to the `Results` tab which displays three major results.

1. Analysis Statistics

Displays a MultiQC consolidated report of overall analysis statistics. This includes FastQC reports and Alignment statistics .

2. Differentially Expressed Genes

On clicking this button a pop up with the list of differentially expressed genes will be displayed.

3. Plots

Displays differentially expressed genes in R plots such as MA plot, Heatmap, PCA plot and box plot.

  

## Walkthrough of pre-constructed pipeline

#### 1. List of pipelines

- WGS data analysis
Enables the user to identify variants from raw sequencing reads and functionally annotate them. Multiple tools are integrated in the WGS analysis pipeline.
- RNA-Seq data analysis
RNA Seq part enables quantification of gene expression. The pipeline provides a final output of a list of differentially expressed genes.
(to be constructed as provided in the CANEapp manual)

  

#### 2. Description of the tools used

Table showing the tools incorporated in iCOMIC

 
|Function  | DNA-Seq Tools |RNA-Seq Tools   |
|--|--|--|
| Quality Control | FastQC, MultiQC, Cutadapt | FastQC, MultiQC, Cutadapt |
|Alignment|GEM-Mapper, BWA-MEM, Bowtie2|Bowtie2, STAR, HISAT22, Salmon|
|Variant Calling|GATK HC, samtools mpileup, FreeBayes, MuSE, GATK Mutect2|--|
| Annotation | Annovar, SnpEff | - |
| Expression Modeller | - | Stringtie, HTSeq,STAR |
| Differential Expression | - | Deseq2, Ballgown |

  - Parameter list for each tool along with tool description
  
- Tools for Quality Control
	1. FastQC
	It is a popular tool that can be used to provide an overview of the basic quality control metrics for raw next generation sequencing data. There are a number different analyses (called modules) that may be performed on a sequence data set. It provides summary graphs enabling the user to decide on the directions for further analysis.
	2.  MultiQC
		MultiQC is a modular tool to aggregate results from bioinformatics analyses across multiple samples into a single report. It collects numerical stats from different modules enabling the user to track the behavior of the data in an efficient manner.
	3.  Cutadapt
  		Cutadapt is a trimming tool that enables the user to remove adapter and primer sequences in an error-tolerant manner. It can also aid in demultiplexing, filtering and modification of single-end and paired-end reads.
- Aligners
	1.  GEM-Mapper
    It is a high-performance mapping tool that performs alignment of sequencing reads against large reference genomes. GEM Mapper has been identified as an efficient mapping tool by a benchmarking analysis performed along with this study. Listed below are some parameters of the tool GEM-Mapper. Others can be found in [GEM-Mapper github page](https://github.com/smarco/gem3-mapper)
    
  |Parameter  | Description |
  |--|--|
  | -t | Threads |
  | -e | --alignment-max-error |
  | --alignment-global-min-identity |Minimum global-alignment identity required |
  | --alignment-global-min-score |Minimum global-alignment score required |

	2.  BWA-MEM
    One of the most commonly used aligners available. It is identified as a faster and accurate algorithm among the algorithms in BWA software package. It is known for aligning long sequence query reads to the reference genome and also performs chimeric alignment. The parameters for BWA-MEM include
	
	|Parameter  | Description |
	|--|--|
	| -t | Threads |
	| -k | minSeedLength |
	| -w |Band width |
	| -d |Z-dropoff |
	| -r |seedSplitRatio |
	|-A |matchScore |
	The other parameters for the tool can be found in [BWA manual page](http://bio-bwa.sourceforge.net/bwa.shtml)
	3.  Bowtie2
    Bowtie2 is a fast and efficient algorithm for aligning reads to a reference sequence. It comprises various modes wherein it supports local, paired-end and gapped alignment. The key parameters for Bowtie2 include:
    
    |Parameter  | Description |
	|--|--|
	| --threads | Threads |
	| --cutoff (n) | Index only the first (n)  bases of the reference sequences (cumulative across sequences) and ignore the rest. |
	| -seed |The seed for pseudo-random number generator. |
	| -N | Sets the number of mismatches to allowed in a seed alignment during [multiseed alignment
	| dvc |the period for the difference-cover sample |
	All parameters for Bowtie2 are listed in [Bowtie2 manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#main-arguments)
	4.  STAR
   STAR is a rapid RNA-Seq read aligner specializing in fusion read and splice junction detection. Important parameters for STAR is given below.
    
    |parameters |Description |
	|----|--|
	|- - runThreadN |NumberOfThreads
	| - -runMode |genomeGenerate
	|- -genomeDir|/path/to/genomeDir|
	|--genomeFastaFiles | /path/to/genome/fasta1 /path/to/genome/fasta2 |
	| --sjdbGTFfile |/path/to/annotations.gtf
	|--sjdbOverhang|ReadLength-1|
	The other parameters for the tool can be found in [STAR manual page](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)

	5.  HISAT2 
	 It is a fast and sensitive alignment program applicable for both RNA seq and Whole-Genome Sequencing data and is known for rapid and accurate alignment of sequence reads to a single reference genome. The key parameters for the tool are given below.
	 
	| parameters |Description |
	|--|--|
	|- x (hisat-idx) |The basename of the index for the reference genome
	| -q |Reads which are FASTQ files
	| --n-ceil (func)|Sets a function governing the maximum number of ambiguous characters (usually Ns and/or .s) allowed in a read as a function of read length.|
	|--ma (int) | Sets the match bonus|
	|--pen-cansplice (int) | Sets the penalty for each pair of canonical splice sites (e.g. GT/AG)| 

	The other parameters for the tool can be found in [HISAT2 manual](http://www.ccb.jhu.edu/software/hisat/manual.shtml)
	6. Salmon
	Acts as a tool for alignment as well as quantifying differential expression. Some of the parameters for the tool include
		 
	| Parameters |Description |
	|--|--|
	|--validateMappings |Enables selective alignment of the sequencing reads when mapping them to the transcriptome|
	|--recoverOrphans|This flag performs orphan “rescue” for reads |
	|--hardFilter|This flag (which should only be used with selective alignment) turns off soft filtering and range-factorized equivalence classes, and removes all but the equally highest scoring mappings from the equivalence class label for each fragment.|
	|--genomeFastaFiles | /path/to/genome/fasta1 /path/to/genome/fasta2 |
	| --numBootstraps |Ables to optionally compute bootstrapped abundance estimates |
	|-p / --threads|The number of threads that will be used for quasi-mapping, quantification, and bootstrapping / posterior sampling (if enabled)|

	The other parameters for the tool can be found in [SALMON Manual Page]		(https://salmon.readthedocs.io/en/latest/salmon.html#mimicbt2)

-  Variant Callers

	1.  GATK HC
	One of the extensively used variant callers. Calls variants from the aligned reds corresponding to the reference genome. Some of the parameters for GATK Haplotype caller are listed beow.
	
	
	| Parameters |Description |
	|----|--|
	|  -contamination|Contamination fraction to filter
	|  -hets|heterozygosity
	|  -mbq|Min base quality score|
	| -minReadsPerAlignStart | Min Reads Per Alignment Start |
	
	
	The complete parameter list is available at [GATK Haplotypecaller article page](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)
	
	2.  Samtools mpileup
	Samtools mpileup together with BCFtools call identifies the variants. Some key parameters to look are listed below. 
	
	| Parameter  | Description |
	|--|--|
	| -d |  --max-depth |
	| -q | Minimum mapping quality for an alignment to be used |
	| -Q |Minimum base quality for a base to be considered |
	Parameters in detail are found in [Samtools-mpileup manual page](http://www.htslib.org/doc/samtools-mpileup.html)
	3.  FreeBayes
	FreeBayes is a variant detector developed to identify SNPs, Indels, MNPs and complex variants with respect to the reference genome. Key parameters for FreeBayes are listed below.
	
	| Parameter  | Description |
	|--|--|
	| -4 | Include duplicate-marked alignments in the analysis. |
	| -m | minimum mapping quality |
	| -q |minimum base quality |
	| -! |minimum coverage |
	| -U |read mismatch limit |
	
	Other parameters can be found in detain in [FreeBayes parameter page](https://vcru.wisc.edu/simonlab/bioinformatics/programs/freebayes/parameters.txt)
	
	4.  GATK Mutect2
	  	This tool identifies somatic mutations in a diseased sample comparing to the provided normal sample. Parameters specific to Mutect2 include 
	  		
	| Parameter  | Description |
	|--|--|
	| --base-quality-score-threshold | Base qualities below this threshold will be reduced to the minimum |
	|--callable-depth | Minimum depth to be considered callable for Mutect stats. Does not affect genotyping. |
	|--max-reads-per-alignment-start |Maximum number of reads to retain per alignment start position. |
	| -mbq |Minimum base quality required to consider a base for calling |
	
	The complete parameter list is available at [GATK Mutect2 manual page](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)
	5.  MuSE
	MuSe is a tool that calls somatic point mutations in normal-tumor sample pairs using a Markov substitution model for evolution.. More information about the tool can be found in[MuSE variant caller tool page](https://bioinformatics.mdanderson.org/public-software/muse/)
- Annotators
	1.  SnpEff
SnpEff tool performs genomic variant annotations and functional effect prediction. Key parameters for the tool SnpEff are listed below.  Detailed list of parameters is given in [SnpEff manual page](http://snpeff.sourceforge.net/SnpEff_manual.html#cmdline)

	| Parameter  | Description |
	|--|--|
	| -t | Use multiple threads|
	| -cancer | perform 'cancer' comparisons (Somatic vs Germline) |
	| -q |Quiet mode|
	| -v |Verbose mode |
	| -csvStats |Create CSV summary file instead of HTML|
	2.  Annovar
Annovar can be used to efficiently annotate functional variants such as SNVs and indels, detected from diverse genomes. The tool also provides the user with multiple annotation strategies namely Gene-based, region-based and filter-based. Key parameters for Annovar include the following.

	| Parameter  | Description |
	|--|--|
	| --splicing_threshold |distance between splicing variants and exon/intron boundary|
	| --maf_threshold| filter 1000G variants with MAF above this threshold |
	| --maxgenethread |max number of threads for gene-based annotation|
	| --batchsize |batch size for processing variants per batch (default: 5m) |
	Details of the tool can be found in [Annovar documentation page](http://annovar.openbioinformatics.org/en/latest/user-guide/gene/)
 -  Expression modellers
	1.  Stringtie
Stringtie is known for efficient and rapid assembly of RNA-Seq alignments into possible transcripts. It employs a novel network flow algorithm and an optional de novo assembly algorithm to assemble the alignments. The important parameters to look into are,

	|  Parameters  |Description |
	|----|--|
	|  --rf |Assumes a stranded library fr-firststrand |
	|  --fr |Assumes a stranded library fr-secondstrand |
	|  --ptf (f_tab)|Loads a list of point-features from a text feature file (f_tab) to guide the transcriptome assembly|
	|  -l (label) |Sets (label) as the prefix for the name of the output transcripts |
	|  -m (int)|Sets the minimum length allowed for the predicted transcripts|
	
	The other parameters for the tool can be found in [STRINGTIE Manual Page](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)

	2.  HTSeq
	It counts the number of mapped reads to each gene. Key parameters for the tool are given below.

	|   Parameters |Description |
	|----|--|
	|  -f |Format of the input data |
	|  -r |For paired-end data, the alignment have to be sorted either by read name or by alignment position |
	|  -s|whether the data is from a strand-specific assay|
	|  -a |skip all reads with alignment quality lower than the given minimum value |
	|  -m|Mode to handle reads overlapping more than one feature|

	The other parameters for the tool can be found in [HTSEQ Manual Page](https://htseq.readthedocs.io/en/release_0.11.1/count.html))
-  Differential Expression tools
	1.  Ballgown
	Assembled transcripts are statistically analysed by the tool. Key arguments for the tool are given below.

	|  Arguments |Description |
	|----|--|
	|  samples |vector of file paths to folders containing sample-specific ballgown data |
	|  dataDir |file path to top-level directory containing sample-specific folders with ballgown data in them|
	|  samplePattern|regular expression identifying the subdirectories of\ dataDir containing data to be loaded into the ballgown object|
	| bamfiles | optional vector of file paths to read alignment files for each sample |
	|  pData |optional data.frame with rows corresponding to samples and columns corresponding to phenotypic variables |
	|  meas|character vector containing either "all" or one or more of: "rcount", "ucount", "mrcount", "cov", "cov_sd", "mcov", "mcov_sd", or "FPKM"|

	The other arguments for the tool can be found in [Ballgown Manual](https://www.bioconductor.org/packages/release/bioc/manuals/ballgown/man/ballgown.pdf)

2.  Deseq2
	Uses negative binomial distribution for testing differential expression. Some of the arguments to look into are given below.
	
	|  Arguments |Description |
	|----|--|
	|  object |A Ranged Summarized Experiment or DESeqDataSet |
	|   groupby |a grouping factor, as long as the columns of object |
	|  run|optional, the names of each unique column in object|
	|  renameCols | whether to rename the columns of the returned object using the levels of the grouping factor|
	|   value |an integer matrix |

	The other arguments for the tool can be found in [DESEQ2 Manual Page](https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf)

## Creating a custom pipeline

#### 1. Shell scripts to be written by the user

#### 2. How to run a custom pipeline

## Viewing and analyzing results

#### 1. MultiQC reports
The tool multiQC compiles the analysis statistics for different tools and provides a consolidated report. The tool is used to visualize the analysis results at multiple stages in iCOMIC. The Quality Control part in iCOMIC analyses the quality of all input reads using FastQC. MultiQC compiles the FastQC reports for each sample and provides a consolidated comprehensive report. MultiQC is used again to summarise the entire analysis statistics. In the case of Whole genome sequencing, the MultiQC report includes a compiled FastQC report, alignment statistics and statistics of the variants identified. On the other hand, the MultiQC report in the RNA Seq analysis part includes results from tools such as FastQC, Cutadapt, and STAR.

#### 2. Plots generated in RNA Seq
Differential Expression tools generate R plots such as MA plot, Heatmap, PCA plot and box plot displays the predicted differentially expressed genes. MA plot helps to find log2 fold changes,Heatmap helps in exploring the count matrix, PCA Plot visualizes the overall effect of experimental covariates and batch effects and box plots used to find count outliers.
#### 3. List of differentially expressed genes in RNA-Seq
#### 4. Variants Called in DNA Seq
iCOMIC displays the variants identified in `vcf` format. In the results tab, the user can click on the `click` button and a pop-up with the `vcf` file will be displayed. 
#### 5. Annotated variants in DNA Seq
Here the `vcf` file of annotated variants are displayed.

## Gallery

#### 1. Screenshots

1. MultiQC Report

2. R plots

3. GUI tabs

#### 2. GIF

1. How to open iCOMIC

2. How to choose a pipeline

3. How to input a file

4. How to choose tools

5. How to choose advanced parameters

6. How to run the pipeline

7. How to view results

## Case examples

#### 1. Installation

#### 2. Configuring iCOMIC

#### 3. Analysis using test data WGS

#### 4. Analysis using test data RNA Seq

#### 5. Results

## Troubleshooting runtime issues

#### 1. FAQs

## Changelog

## Glossary
1.  GTF- Gene Transfer Format
    
2.  NGS - Next Generation Sequencing
    
3.  GUI - Graphical User Interface
    
4.  DNA-Seq - DNA Sequencing (Whole Genome/exome Sequencing)
    
5.  RNA-Seq - RNA Sequencing
    
6.  iCOMIC - Integrating the Context Of Mutations In Cancer 

