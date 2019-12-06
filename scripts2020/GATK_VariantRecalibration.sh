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

## variants-sets
# True sites training resource: HapMap
truetraining15=reference/hg38_v0_hapmap_3.3.hg38.vcf.gz

# True sites training resource: Omni
truetraining12=reference/hg38_v0_1000G_omni2.5.hg38.vcf.gz

# Non-true sites training resource: 1000G
nontruetraining10=reference/hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# Known sites resource, not used in training: dbSNP
knowntraining2=reference/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz

# indels True sites training resource: Mills
truetrainingindel12=reference/hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# indels Axiom
axiom10=reference/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz

# path to GATK executable
GATK=/opt/biotools/gatk

#raw_variants_vcf=variants/NA12878.vcf.gz
raw_variants_vcf=variants/NA12878-chr22.vcf.gz

prefix=variants/$(basename ${raw_variants_vcf%.vcf.gz})

##############
##   VQSR   ##
##############

# Filter the VCF file according to excess heterogeneity
# Threshold for marking calls as ExcessHet 
# if the probability that the number of heterozygotes is higher than conform with the law of Hardy and Weinberg, 
# expressed as - 10 log (p). The default value of 54.69 corresponds to a p-value of 3.4 10-6 and a z-zcore of -4.5.
# This produces a VCF callset where any record with ExcessHet greater than 54.69 is filtered with the ExcessHet label in the FILTER column.

threshold=54.69
java -Xmx24g -Xms24g -jar $GATK/gatk.jar \
    VariantFiltration \
    --variant ${raw_variants_vcf} \
    --filter-expression "ExcessHet > ${threshold}" \
    --filter-name "ExcessHet" \
    -O ${prefix}_excesshet_filtered.vcf

# extract a 6-column version of the data for recalibration
java -Xmx8G -jar $GATK/gatk.jar \
    MakeSitesOnlyVcf \
    -I ${prefix}_excesshet_filtered.vcf \
    -O ${prefix}_excesshet_sitesonly.vcf

# Build the SNP recalibration model
java -Xmx24g -Xms24g -jar $GATK/gatk.jar \
    VariantRecalibrator \
    -R ${reference_fa} \
    -V ${prefix}_excesshet_sitesonly.vcf \
    -O ${prefix}_recalibrate_SNP.recal.vcf \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${truetraining15} \
    --resource:omni,known=false,training=true,truth=true,prior=12.0 ${truetraining12} \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${nontruetraining10} \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${knowntraining2} \
    --use-annotation DP \
    --use-annotation QD \
    --use-annotation FS \
    --use-annotation SOR \
    --use-annotation MQ \
    --use-annotation MQRankSum \
    --use-annotation ReadPosRankSum \
    --mode SNP \
    --tranches-file ${prefix}_recalibrate_SNP.tranches \
    --tranche 100.0 \
    --tranche 99.9 \
    --tranche 99.0 \
    --tranche 90.0 \
    --rscript-file ${prefix}_recalibrate_SNP_plots.R

# Build the Indel recalibration model
java -Xmx24g -Xms24g -jar $GATK/gatk.jar \
    VariantRecalibrator \
    -R ${reference_fa} \
    -V ${prefix}_excesshet_sitesonly.vcf \
    -O ${prefix}_recalibrate_indels.recal.vcf \
    --resource:mills,known=false,training=true,truth=true,prior=12 ${truetrainingindel12} \
    --resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiom10} \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2 ${knowntraining2} \
    --trust-all-polymorphic \
    --use-annotation QD \
    --use-annotation DP \
    --use-annotation FS \
    --use-annotation SOR \
    --use-annotation MQRankSum \
    --use-annotation ReadPosRankSum \
    --mode INDEL \
    --max-gaussians 4 \
    --tranches-file ${prefix}_recalibrate_indels.tranches \
    --tranche 100.0 \
    --tranche 99.95 \
    --tranche 99.9 \
    --tranche 99.5 \
    --tranche 99.0 \
    --tranche 97.0 \
    --tranche 96.0 \
    --tranche 95.0 \
    --tranche 94.0 \
    --tranche 93.5 \
    --tranche 93.0 \
    --tranche 92.0 \
    --tranche 91.0 \
    --tranche 90.0 \
    --rscript-file ${prefix}_recalibrate_INDEL_plots.R

# Apply the desired level of recalibration to the SNPs in the call set
java -Xmx24g -Xms24g -jar $GATK/gatk.jar \
    ApplyVQSR \
    -R ${reference_fa} \
    -V ${prefix}_excesshet_filtered.vcf \
    --mode SNP \
    --truth-sensitivity-filter-level 99.0 \
    --recal-file ${prefix}_recalibrate_SNP.recal.vcf \
    --tranches-file ${prefix}_recalibrate_SNP.tranches \
    -O ${prefix}_recalibrated_snps_raw_indels.vcf.gz

# Apply the desired level of recalibration to the Indels in the call set
java -Xmx24g -Xms24g -jar $GATK/gatk.jar \
    ApplyVQSR \
    -R ${reference_fa} \
    -V ${prefix}_recalibrated_snps_raw_indels.vcf.gz \
    --recal-file ${prefix}_recalibrate_indels.recal.vcf \
    --tranches-file ${prefix}_recalibrate_indels.tranches \
    --truth-sensitivity-filter-level 99.7 \
    --create-output-variant-index true \
    -mode INDEL \
    -O ${prefix}_VQSR.vcf.gz

