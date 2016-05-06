#! /usr/bin/env bash
# script varscan2-mpileup2cns_call-variants.sh

# required
# samtools version: 1.x (htslib)
# varscan2 version 2.4.1
# vcftools 0.1.14
# vcf2index custom function

# make sure you have defined VARSCAN in your .bashrc
# export VARSCAN=<edit here the path to VarScan.v2.4.1.jar>

# if not done yet, source user-defined .myfunctions (including vcf2index)
# . $HOME/.myfunctions
# vcf2index is reproduced here for training
vcf2index {
# takes a vcf file as argument ($1)
# keeps header and sorts remaining lines
# compresses sorted file with bgzip
# indexes it with tabix for vcftools applications
# author: Stephane Plaisance
# version: 1.0
if [ $# -eq 1 ]; then
	title=$(basename "$1" ".vcf")
	( grep ^"#" $1 | perl -pi -e "s/$title.*$/$title/";
	grep -v ^"#" $1 | sort -k 1V,1 -k 2n,2 ) | bgzip -c > $1".gz"
	tabix -p vcf $1".gz"
else
	echo "usage: vcf2index <vcf-file>"
fi
}

# we call from the full mapping results for chr21
# calling from the low coverage training data would not give confident calls
infolder=hg19_bwa-mapping
infile=shuffled_PE_NA18507_GAIIx_100_chr21_mem_mdup.bam

outfolder=varscan2_variants
mkdir -p ${outfolder}

ref=ref/HiSeq_UCSC_hg19.fa

# call SNV and InDels variants from samtools mpileup

## limit allocated RAM
maxmem="4G"

samtools mpileup -f ${ref} ${infolder}/${infile} | \
	java -Xmx${maxmem} -jar $VARSCAN/VarScan.v2.4.1.jar mpileup2cns \
	--variants --output-vcf 1 --p-value 0.01 \
	>${outfolder}/NA18507_mpileup2cns.vcf \
	2>${outfolder}/varscan2_mpileup2cns.err

## compress and index result
bgzip -c ${outfolder}/NA18507_mpileup2cns.vcf \
	> ${outfolder}/NA18507_mpileup2cns.vcf.gz &&
	tabix -p vcf ${outfolder}/NA18507_mpileup2cns.vcf.gz

# extract the chr21 subset of calls
vcftools --gzvcf ${outfolder}/NA18507_mpileup2cns.vcf.gz \
	--chr chr21 \
	--out ${outfolder}/chr21_NA18507_varscan \
	--recode

# recoding introduces the word recode at the end of the file-name
# remove the recode tag from the name
mv ${outfolder}/chr21_NA18507_varscan.recode.vcf \
	${outfolder}/chr21_NA18507_varscan.vcf

# sort, compress and index using our custom 'vcf2index' function)
vcf2index ${outfolder}/chr21_NA18507_varscan.vcf && \
	rm ${outfolder}/chr21_NA18507_varscan.vcf

echo