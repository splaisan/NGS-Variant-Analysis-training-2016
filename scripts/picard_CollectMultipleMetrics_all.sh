#! /usr/bin/env bash
# script: 'picard_CollectMultipleMetrics.sh'
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


# OPT: PROGRAM={CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, 
#                       QualityScoreDistribution, MeanQualityByCycle} 


echo "# processing ${infile}" 
# run in parallel if the computer has enough capacity
( java -Xmx4g -jar $PICARD/picard.jar CollectMultipleMetrics \
	R=${ref} \
	I=${infolder}/${infile} \
	O=${outfolder}/${infile%%.bam} \
	PROGRAM=CollectAlignmentSummaryMetrics \
	PROGRAM=CollectInsertSizeMetrics \
	PROGRAM=QualityScoreDistribution \
	PROGRAM=MeanQualityByCycle \
	2>${outfolder}/${infile%%.bam}_CollectMultipleMetrics.err )

echo

