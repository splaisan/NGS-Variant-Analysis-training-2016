#! /usr/bin/env bash
## script: 'picard_CollectGcBiasMetrics.sh'
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

# send errors log file with '2>'
java -jar $PICARD/picard.jar CollectGcBiasMetrics \
	VERBOSITY=ERROR \
	R=${ref} \
	I=${infolder}/${infile} \
	O=${outfolder}/${infile%%.bam}_GCbias.txt \
	CHART=${outfolder}/${infile%%.bam}_GCbias.pdf \
	S=${outfolder}/${infile%%.bam}_GCbias_summary.txt \
	2>${outfolder}/${infile%%.bam}_CollectGcBiasMetrics.err
echo

