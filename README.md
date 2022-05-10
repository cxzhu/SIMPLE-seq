# SIMPLE-seq

Processing of SIMPLE-seq datasets.

#### Please have the following softwares installed first:

- bowtie, http://bowtie-bio.sourceforge.net/index.shtml
   
   For bowtie version=1.x, please modify <code>shellscrips/01.pre_process_paired_tag_fastq.sh</code> according to the comments in that file.

- bowtie2, http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

- samtools, http://www.htslib.org/
   samtools version >= 1.3.1 is required.

- Trim_galore, https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

- Optional: FastQC, https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

- Compile the "simplecov" tool:

	 <code>cd cimpleconv</code>

	 <code>sh make.sh</code>



