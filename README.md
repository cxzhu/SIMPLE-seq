# SIMPLE-seq

Processing of SIMPLE-seq datasets.

SIMPLE-seq is a scalable method for joint analysis of 5mC and 5hmC from single cells. This repository provide the scripts for decoding the cellular barcodes of SIMPLE-seq datasets (modified from ligation-based combinatorial barcoding from SPLiT-seq), and for the identification of 5mC and 5hmC sites for individual cells.

#### Please have the following softwares installed first:

- bowtie, http://bowtie-bio.sourceforge.net/index.shtml
   
- bowtie2, http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

- samtools, http://www.htslib.org/
   samtools version >= 1.3.1 is required.

- Trim_galore, https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

- Optional: FastQC, https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

- Compile the "simplecov" tool:

	 <code>cd simpleconv</code>

	 <code>sh make.sh</code>


#### Analysis of SIMPLE-seq datasets include the following steps:


## 1. Pre-processing
Extract cellular barcode from Read2, map the reads to reference cell_ID, and convert the mapped cell ID samfiles to useable fastq files.

*** Please modification the paths to reference files according to the annotations in the script file.

Use <code>sh shellscrips/01.pre_process_simple_seq_fastq.sh [sample_prefix]</code>.

#### The output of this step includes:

1. <code>Sample_combined.fq.gz</code>  This file is a combined fastq file including Read1 sequences/qualities and barcode sequences extracted from Read2.

2. <code>Sample_BC.sam</code> This is a temporally file used to assign extracted barcode sequences to Cellular Barcode. Please delete this file if you have successful obtained <code>Sample_BC_cov.fq.gz</code>.

3. <code>Sample_BC_cov.fq.gz</code> This is the fastq file with Read1 sequences and qualities, the Cellular Barocde and UMI from Read2 are now in ReadName section of the fastq file (and subsequent alignment files).

## 2. Mapping to the genome
As SIMPLE-seq only introduce "C-to-T" mutations on 5mC and 5hmC sites, we used bowtie2 (instead of other methylation aligner) for mapping.

Use <code>sh shellscrips/02.proc_mapping.sh [sample_prefix]</code>.


## 3. Split the alignment files to 5mC and 5hmC
This step is to split 5mC and 5hmC reads to seperate alignment files (bam files) based on the indicator sequences.

Use <code>perl perlscripts/02.split_modality.pl [sample_sorted.bam]</code>.

Three files will be generated, including <code>[sample_sorted.bam_5mC.bam], [sample_sorted.bam_5hmC.bam] and [sample_sorted.bam_other.bam] </code>. Reads cannot be perfectly assigned to 5mC or 5hmC will be written to "XXX_other.bam".

## 4. Generate cell-to-modification abundance matrices
This step will convert the bam files to an intermediate modification information file and then generate cell-to-modification abundance matrices.

Step.1 <code>perl perlscripts/03.bam2srf.pl [sample_sorted.bam_5mC/5hmC.bam]</code>.

Step.2 <code>perl perlscripts/04.srf2mtx.pl [input.rsf] [binsize]</code>.

The resulting matrix can be used for downstream single-cell analysis.
