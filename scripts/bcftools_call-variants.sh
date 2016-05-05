#! /usr/bin/env bash
# script bcftools_call-variants.sh

# required
# samtools v1.x (htslib)
# bcftools v1.x (htslib)
# vcftools 0.1.14
# vcf2index custom function

# if not done yet, source user-defined functions (including vcf2index)
# . /home/bits/.myfunctions

function vcf2index() {
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
name=${infile%%.bam}
outfolder=bcftools_htslib_variants
mkdir -p ${outfolder}

ref=ref/HiSeq_UCSC_hg19.fa

# call variants from pileup, filter and convert calls to VCF format
# samtools mpileup arguments
#  -u generate uncompress BCF output (to feed the pipe)
#  -f reference genome in fasta format related to the BWA index
# bcftools call arguments
#  -v output potential variant sites only (force -c)
#  -c SNP calling (force -e)
#  -O u for bcf-uncompressed output
# exclude calls where more than 1000 reads support a variation

## REM: the next command is made of two parts contained between '(...)' with separate error logs (2>...)
## only if (part-1) succeeds ('&&') => then (part-2) is started.

## part-1 creates a raw BCF file
(samtools mpileup -uf ${ref} ${infolder}/${infile} | \
	bcftools call -vc -O u > ${outfolder}/NA18507_var_bcftools.raw.bcf \
	2>${outfolder}/samtools_mpileup_bcftools_raw.err) && \

## part-2 processes the BCF data, filters it, and creates a VCF dataset
(bcftools view ${outfolder}/NA18507_var_bcftools.raw.bcf | \
	vcfutils.pl varFilter -D1000 > \
	${outfolder}/NA18507_var_bcftools.flt-D1000.vcf \
	2>${outfolder}/samtools_mpileup_bcftools_filter.err)

# Sort, compress, and index the variant VCF files
vcf2index ${outfolder}/NA18507_var_bcftools.flt-D1000.vcf && \
	rm ${outfolder}/NA18507_var_bcftools.flt-D1000.vcf

# extract the chr21 subset of calls
vcftools --gzvcf ${outfolder}/NA18507_var_bcftools.flt-D1000.vcf.gz \
	--chr chr21 \
	--out ${outfolder}/chr21_NA18507_var_bcftools.flt-D1000 \
	--recode

# recoding introduces the word recode at the end of the file-name
# remove the recode tag from the name to simplify
mv ${outfolder}/chr21_NA18507_var_bcftools.flt-D1000.recode.vcf \
	${outfolder}/chr21_NA18507_var_bcftools.flt-D1000.vcf

# sort, compress and index using our custom 'vcf2index' function)
vcf2index ${outfolder}/chr21_NA18507_var_bcftools.flt-D1000.vcf && \
	rm ${outfolder}/chr21_NA18507_var_bcftools.flt-D1000.vcf

echo