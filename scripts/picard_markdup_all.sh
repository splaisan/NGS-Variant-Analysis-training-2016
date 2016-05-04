#! /usr/bin/env bash
## script: 'picard_markdup.sh'
## ©SP-BITS, 2013 v1.0
# last edit: 2016-02-19
 
# required: 
# Picard Version: 2.1.0

# 10% sample
# infile="hg19_bwa-mapping/shuffled_10pc_PE_NA18507_GAIIx_100_chr21_mem-pe.bam"
 
# full data
infile="hg19_bwa-mapping/shuffled_PE_NA18507_GAIIx_100_chr21_mem-pe.bam"
 
echo "# marking duplicates in ${infile}"

# we extract part of the infile and replace the end of it by a new text and extension
# ${infile%%-pe.bam} removes the text '-pe.bam' from the end of 'infile'

java -jar $PICARD/picard.jar MarkDuplicates \
        I=${infile} \
        O=${infile%%-pe.bam}_mdup.bam \
        ASSUME_SORTED=TRUE \
        REMOVE_DUPLICATES=FALSE \
        CREATE_INDEX=TRUE \
        M=${infile%%-pe.bam}_duplicate_metrics.txt \
        VALIDATION_STRINGENCY=LENIENT \
        2>${infile%%-pe.bam}_MarkDuplicates.err
echo
 
echo "# results are"
head -8 ${infile%%-pe.bam}_duplicate_metrics.txt| tail -2 | transpose -t | column -t
echo
