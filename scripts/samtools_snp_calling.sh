#!/usr/bin/env bash

module load apps/gcc/5.2/samtools

REF_FILE=data/ref_files/Parus_major_1.04.chrLGE22.fa
VCF_OUT=vcf_files

if [ ! -d "$VCF_OUT" ]; then
    mkdir -p ${VCF_OUT} # create the directory where the vcf files are to be written
fi

ls data/bam_files/*.bam > bams.list # create input list of bam files for samtools mpileup

bcftools_1_3=/usr/local/extras/Genomics/apps/bcftools/1.3/bin/bcftools

samtools mpileup -ugf ${REF_FILE} -b bams.list | ${bcftools_1_3} call -vmO z -f GQ -o ${VCF_OUT}/samtools.chrLGE22.raw.snps.indels.vcf.gz

