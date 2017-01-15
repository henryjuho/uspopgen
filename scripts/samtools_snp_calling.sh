#!/usr/bin/env bash

module load apps/gcc/5.2/samtools

REF_FILE=data/ref_files/Parus_major_1.04.chrLGE22.fa
VCF_OUT=vcf_files

if [ ! -d "$VCF_OUT" ]; then
    mkdir -p ${VCF_OUT} # create the directory where the vcf files are to be written
fi

ls data/bam_files/*.bam > bams.list # create input list of bam files for samtools mpileup

export PATH=/usr/local/extras/Genomics/apps/bcftools/1.3/bin/:$PATH

samtools mpileup -q 20 -Q 10 -ugf  ${REF_FILE} -b bams.list | bcftools call -vmO z -f GQ -o ${VCF_OUT}/samtools.chrLGE22.raw.snps.indels.vcf.gz

bcftools index -t ${VCF_OUT}/samtools.chrLGE22.raw.snps.indels.vcf.gz 
