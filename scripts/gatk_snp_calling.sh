#!/usr/bin/env bash

module load apps/binapps/GATK

REF_FILE=data/ref_files/Parus_major_1.04.chrLGE22.fa
VCF_OUT=vcf_files

if [ ! -d "$VCF_OUT" ]; then
    mkdir -p ${VCF_OUT} # create the directory where the vcf files are to be written
fi

ls data/gvcf_files/*.vcf > gvcfs.list # create input list of gvcf files for GenotypeGVCFs

java -Xmx3g -jar $GATKHOME/GenomeAnalysisTK.jar \
     -T GenotypeGVCFs \
     -R ${REF_FILE} \
     -V gvcfs.list \
     -o ${VCF_OUT}/gatk.chrLGE22.raw.snps.indels.vcf.gz
