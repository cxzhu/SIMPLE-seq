s=$1 ## sample prefix
ref="path-to-your-bowtie2-reference-genome" ## modify this to your reference directory

trim_galore ${s}_cov.fq.gz

bowtie2 -x $ref -U ${s}_cov_trimmed.fq.gz -p 4 | gzip - > ${s}_mapped.sam.gz

zcat ${s}_mapped.sam.gz | samtools sort - -o ${s}_mapped_sorted.bam
