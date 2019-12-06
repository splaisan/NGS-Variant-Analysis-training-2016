#!/bin/env bash

# from https://software.broadinstitute.org/gatk/documentation/article.php?id=2805
# https://gatkforums.broadinstitute.org/gatk/discussion/23216/how-to-filter-variants-either-with-vqsr-or-by-hard-filtering
# SP@VIB-NC, 2019-12-06
# compatible with GATK 4.1.4.1

workdir=/data/NC_projects/BITS_variants_2020
cd ${workdir}

# recalibration sources
## genome
reference_fa=reference/Homo_sapiens.UCSC.hg38.fa

# path to GATK executable
GATK=/opt/biotools/gatk

# set below to match your RAM availability
javaopts="-Xmx24g -Xms24g"

raw_variants_vcf=variants/NA12878-chr22.vcf.gz

prefix=variants/$(basename ${raw_variants_vcf%.vcf.gz})

########################
##   HARD filtering   ##
########################

###########################################
# 1) hard-Filter excess heterogeneity
###########################################

# Threshold for marking calls as ExcessHet 
# if the probability that the number of heterozygotes is higher than conform with the law of Hardy and Weinberg, 
# expressed as - 10 log (p). The default value of 54.69 corresponds to a p-value of 3.4 10-6 and a z-zcore of -4.5.
# This produces a VCF callset where any record with ExcessHet greater than 54.69 is filtered with the ExcessHet label in the FILTER column.

threshold=54.69
java ${javaopts} -jar $GATK/gatk.jar \
    VariantFiltration \
    --variant ${raw_variants_vcf} \
    --filter-expression "ExcessHet > ${threshold}" \
    --filter-name "ExcessHet" \
    -O ${prefix}_excesshet_filtered.vcf.gz

###########################################
# 2) hard-Filter SNPs on multiple metrics
###########################################

# produces a VCF with records with SNP-type variants only.
java ${javaopts} -jar $GATK/gatk.jar \
    SelectVariants \
    -V ${raw_variants_vcf} \
    --select-type-to-include SNP \
    -O ${prefix}_snp.vcf.gz

# This produces a VCF with the same variant records now annotated with filter status. Specifically, if a record passes all the filters, it receives a PASS label in the FILTER column. 
# A record that fails a filter #receives the filter name in the FILTER column, e.g. SOR3. 
# If a record fails multiple filters, then each failing filter name appears in the FILTER column separated by semi-colons ; e.g. "MQRankSum-12.5;ReadPosRankSum-8".

java ${javaopts} -jar $GATK/gatk.jar \
    VariantFiltration \
    -V ${prefix}_snp.vcf.gz \
    --filter "QD < 2.0" --filter-name "QD2" \
    --filter "QUAL < 30.0" --filter-name "QUAL30" \
    --filter "SOR > 3.0" --filter-name "SOR3" \
    --filter "FS > 60.0" --filter-name "FS60" \
    --filter "MQ < 40.0" --filter-name "MQ40" \
    --filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${prefix}_snp_filtered.vcf.gz

###########################################
# 3) hard-Filter INDELs on multiple metrics
###########################################

# produces a VCF with records with INDEL-type variants only.
java ${javaopts} -jar $GATK/gatk.jar \
    SelectVariants \
    -V ${raw_variants_vcf} \
    --select-type-to-include INDEL \
    -O ${prefix}_indels.vcf.gz

java ${javaopts} -jar $GATK/gatk.jar \
    VariantFiltration \
    -V ${prefix}_indels.vcf.gz \
    --filter "QD < 2.0" --filter-name "QD2" \
    --filter "QUAL < 30.0" --filter-name "QUAL30" \
    --filter "FS > 200.0" --filter-name "FS200" \
    --filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O ${prefix}_indels_filtered.vcf.gz

###########################################
## 4) merge SNP and Indel filtered calls    
###########################################

# conbine filtered into one file using Picard
java ${javaopts} -jar $GATK/gatk.jar \
    MergeVcfs \
    -I ${prefix}_snp_filtered.vcf.gz \
    -I ${prefix}_indels.vcf.gz \
    -R ${reference_fa} \
    -O ${prefix}_snp_indel_filtered.vcf.gz
    
