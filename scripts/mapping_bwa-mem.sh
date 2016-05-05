#! /usr/bin/env bash
## script: 'mapping_bwa-mem.sh'
## ©SP-BITS, 2013 v1.1
# last edit: 2016-02-18

# required:
# bwa version: 0.7.12-r1039
# Picard Version: 2.1.0
# bamutils Version: 1.0.12b

## define global variables
refgen=ref/HiSeq_UCSC_hg19.fa

## select one of the folowing three blocks and comment the other two
infolder=reads

## full data
# f=shuffled_NA18507_GAIIx_100_chr21
# pfx=all

## 10% sample
f=shuffled_10pc_PE_NA18507_GAIIx_100_chr21
pfx=sample

# create new folder
outfolder=hg19_bwa-mapping
mkdir -p ${outfolder}

# reads
fq1=${infolder}/${f}_2_1.fq.gz
echo "# forward reads: ${fq1}"
fq2=${infolder}/${f}_2_2.fq.gz
echo "# reverse reads: ${fq2}"

# marking secondary hits with -M to comply with Picard
# using 'nthr' processors in parallel (again limited by our RAM!)
# mem is more demanding and needs more than 3Gb per thread
nthr=1
echo "# using ${nthr} processor threads"

rgstring='@RG\tID:NA18507\tLB:lib-NA18507\tPU:unkn-0.0\tPL:ILLUMINA\tSM:GAIIx-chr21-BWA.mem'

echo
echo "# mapping paired reads with **bwa mem**"

## default numeric settings are left unchanged as:
#  -k INT     minimum seed length [19]
#  -w INT     band width for banded alignment [100]
#  -d INT     off-diagonal X-dropoff [100]
#  -r FLOAT   look for internal seeds inside a seed longer than {-k}
#             * FLOAT [1.5]
#  -c INT     skip seeds with more than INT occurrences [10000]
#  -A INT     score for a sequence match [1]
#  -B INT     penalty for a mismatch [4]
#  -O INT     gap open penalty [6]
#  -E INT     gap extension penalty; a gap of size k cost {-O} + {-E}*k [1]
#  -L INT     penalty for clipping [5]
#  -U INT     penalty for an unpaired read pair [17]
#  -T INT     minimum score to output [30]
# '-p': first query file consists of interleaved paired-end sequences
#
# '-M': mark shorter split hits as secondary (for Picard/GATK compatibility)

bwape=${f}_mem-pe

# store the full command line to include it in the next part
cmd="bwa mem -R \"${rgstring}\" \
	-M \
	-t ${nthr} \
	${refgen} \
	${fq1} ${fq2} > ${outfolder}/${bwape}.sam"

# show and execute the command
echo "# ${cmd}"
eval ${cmd} 2>${outfolder}/${pfx}_bwa_mem.err

########################### post process results ##############
#### add @PG group to results with **bamUtil 'bam polishBam' **
# important formatting issues here
# $'\t' add a 'tab' character
# $'\'' add single quotes around ${cmd} to avoid it being interpreted
# finally, $pgstring is called quoted to be replaced by its building variables

# we first detect the current version of BWA to store it
bwaver=$(expr "$(bwa 2>&1)" : '.*Version:\ \(.*\)Contact.*')

# we then create a '@PG' line to store our action and include the 'BWA version number' AND the 'full command'
# clean extra spaces with sed and the POSIX character class [:space:]
# some explanation: [:space:] is functionally identical to [ \t\r\n\v\f]
cmd=$(echo ${cmd} | sed -e "s/[[:space:]]\+/ /g")

pgstring=@PG$'\t'ID:BWA$'\t'PN:bwa$'\t'VN:${bwaver}$'\t'CL:$'\''${cmd}$'\''

# clean up initial file after success to save space
bam polishBam --PG "${pgstring}" \
	-i ${outfolder}/${bwape}.sam \
	-o ${outfolder}/${bwape}_pg.sam \
	-l ${outfolder}/polishBam_${bwape}.log && rm ${outfolder}/${bwape}.sam

##### convert to sorted bam, index & cleanup
java -jar $PICARD/picard.jar SortSam \
	I=${outfolder}/${bwape}_pg.sam \
	O=${outfolder}/${bwape}.bam \
	SO=coordinate \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=LENIENT && rm ${outfolder}/${bwape}_pg.sam

##### get basic stats from the resulting bam file
samtools flagstat ${outfolder}/${bwape}.bam \
	>${outfolder}/${bwape}_flagstat.txt && \
	cat ${outfolder}/${bwape}_flagstat.txt

echo "# finished mapping ${f} reads with BWA mem"

