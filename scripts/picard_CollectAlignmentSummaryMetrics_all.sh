#! /usr/bin/env bash
## script: 'picard_CollectAlignmentSummaryMetrics.sh'
## ©SP-BITS, 2013 v1.0
# last edit: 2016-02-19
# Picard Version: 2.1.0
 
infolder=hg19_bwa-mapping
outfolder=bwa-mappingQC
mkdir -p ${outfolder}
 
ref=ref/HiSeq_UCSC_hg19.fa
 
# 10% sample
# infile=shuffled_10pc_PE_NA18507_GAIIx_100_chr21_mem_mdup.bam
 
# full data
infile=shuffled_PE_NA18507_GAIIx_100_chr21_mem_mdup.bam
 
echo "# processing ${infile}"
 
# collect information
(java -jar $PICARD/picard.jar CollectAlignmentSummaryMetrics \
	R=${ref} \
	I=${infolder}/${infile} \
	O=${outfolder}/${infile%%.bam}_AlignmentSummaryMetrics.txt \
	2>${outfolder}/${infile%%.bam}_CollectAlignmentSummaryMetrics.err)
 
# ignore header and empty lines and rotate the remaining content
grep -v -E '(^#|^$)' \
	${outfolder}/${infile%%.bam}_AlignmentSummaryMetrics.txt | \
	transpose -t | \
	column -t > \
	${outfolder}/tr_${infile%%.bam}_AlignmentSummaryMetrics.txt
echo
