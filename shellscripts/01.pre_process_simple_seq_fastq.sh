s=$1 # with ${s}_R1.fq.gz and ${s}_R2.fq.gz in the working directory



## Step 1
## (a) for PBMC and mESC with 19-bp barcode
p="path-to-cell-id-reference"  # reference file is in "Cell_BC_reference/" of github directory, "cell_id" for mESC and PBMC,
simpleconv combine ${s} 

## (b) for brain with 21-bp barcode
## p="path-to-cell-id-reference"  # reference file is in "Cell_BC_reference/" of github directory,  "SIMPLE_Brain_ID" for brain.
## simpleconv combine3 ${s} 

### version for bowtie will chage for use of this step
zcat ${s}_combined.fq.gz | bowtie ${p} - --norc -m 1 -v 1 -S ${s}_BC.sam ## for bowtie 0.x
# zcat ${s}_combined.fq.gz | bowtie -x ${p} - --norc -m 1 -v 1 -S ${s}_BC.sam ## for bowtie 1.x

#### This step convert to Celluar Barcode mapped reads to fastq files.
#### may need to mofify this script according to read length format (GEO/Illumina/BGI)
perl perlscripts/01.sam2fastq.pl ${s}

