#! /usr/bin/env bash
## script: 'run_fastqc.sh'
## ©SP-BITS, 2013 v1.1
# last edit: 2016-02-18
 
# required: 
# FastQC

infolder="reads"

# full read files
inprefix="shuffled_PE_NA18507_GAIIx_100_chr21"

# 10%-sample read files
# inprefix="shuffled_10pc_PE_NA18507_GAIIx_100_chr21"

# create folder for QC-results inside the read folder
outfolder=readQC
mkdir -p ${infolder}/${outfolder}

# on each reads file
for fq in ${infolder}/${inprefix}*.fq.gz; do
	echo "analysing ${fq}"
	# perform QC test using 1 thread and zip results
	fastqc -t 1 -o ${infolder}/${outfolder} --noextract ${fq}

done