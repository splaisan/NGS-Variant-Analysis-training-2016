#! /usr/bin/env bash
## script: 'qualimap_bamqc.sh'
## ©SP-BITS, 2014 v1.0
# last edit: 2014-04-09

# required:
# qualimap 0.8
# R dependencies

# process all bam files in 'hg19_bwa-mapping'
infolder=hg19_bwa-mapping
outfolder=bwa-mapping-qualimapQC
mkdir -p ${outfolder}

mygtf=ref/Homo_sapiens.hg19_chr21.gtf

for b in ${infolder}/*_mdup.bam; do
outname=$(basename ${b} .bam)

echo "# processing ${b}"

# read the full help to find more about other parameters
# we output here to HTML but PDF is also available for portability and compactness
# we restrict the analysis to chr21 using annotation data from the Qualimap site
## Human genome annotations from Ensembl database (v.64)
qualimap bamqc \
	-bam ${b} \
	-gd HUMAN \
	-c \
	-gff ${mygtf} \
	-outdir ${outfolder}/${outname}_qualimap-results \
	-outformat HTML

# for PDF, replace last line by
#	-outformat PDF
#	-outfile ${outname}_chr21_report.pdf \

done