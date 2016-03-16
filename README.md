# SNP calling and Filtering

We will use some whole genome data from the great tit to look at SNP calling and filtering.
The data we will 10 great tit individuals from Europe, they are a subset of the 29 birds that were sequenced and 
in [Laine et al. (2016)](http://www.nature.com/ncomms/2016/160125/ncomms10474/full/ncomms10474.html).

The 10 birds were sequenced with illumina paired-end data with 100bp reads to a mean depth of between 9-14x.
In the Laine et al. (2016) paper the SNP calling was performed using GATK, Platypus and ANGSD.

## Programs required
All the programs we will use are available on iceberg through the module system or from the genomics repository.

The following programs will be used:
    
* samtools v1.2 (available through the module system)
* bcftools v1.3 (available in the genomics repository)
* GATK v3.4 (available through the  module system)
* ANGSD v0.911

## SNP callers

We will call SNPs using GATK and samtools/bcftools.


### GATK SNP calling

    bash scripts/gatk_snp_calling.sh 

This compressed vcf file will be written in the vcf_files folder.

### SAMTOOL/BCFTOOLS SNP calling


    bash scripts/samtools_snp_calling.sh
    
### Understanding the VCF format
First, we will take a look at the VCF format, as this is the output format for most SNP calling and genotyping programs
 
The VCF format is composed of a header section where each line begins with '##' and the headers describing the columns
are located on the line starting with '#chrom'. The columns of the VCF from left to right are the chromosome/scaffold,  
position within the reference chromosome/scaffold, an ID (here always ‘.’), the reference allele, the variant allele,
the SNP/INDEL quality score, the filter field, info field, the sample format field and the samples make up the remaining
columns (10 in this case).

The INFO field containing a lot of annotations for the variant site across samples (e.g Depth, Mapping Qualtity etc.) and
these values may be used for filtering (see below). The format for the genotype information is explained in the header
of VCF file. Each row following the header section is a variant site in the VCF, either a SNP or an indel.

Now, take a look at the GATK file.

    less -S vcf_files/gatk.chrLGE22.raw.snps.indels.vcf.gz 

Note that different SNP callers will have some differences in the annotations present in the INFO field and differences
in the format fields. What differences do you see between the samtools VCF and the GATK VCF files?

### Comparing the output from the two callers

We will count and compare the unfiltered SNP calls using some common command line tools and the bedtools program.

Count the number of variants in the vcf file using zgrep (grep for compressed files).

    zgrep -cv ^# vcf_files/gatk.chrLGE22.raw.snps.indels.vcf.gz 
    zgrep -cv ^# vcf_files/samtools.chrLGE22.raw.snps.indels.vcf.gz 
    
Which tool calls the most variants?

Now we will extract only biallelic SNPs and discard the INDELs from the VCF using GATK's SelectVariants tool
    
    module load apps/binapps/GATK
    java  -Xmx3g -jar $GATKHOME/GenomeAnalysisTK.jar -T SelectVariants -R data/ref_files/Parus_major_1.04.chrLGE22.fa -V vcf_files/gatk.chrLGE22.raw.snps.indels.vcf.gz -o vcf_files/gatk.chrLGE22.raw.snps.vcf.gz -selectType SNP -restrictAllelesTo BIALLELIC 

We will extract the SNPs from the samtools vcf, but first we need to index it, as the GATK tools require this.

    tabix -p vcf vcf_files/samtools.chrLGE22.raw.snps.indels.vcf.gz
    java  -Xmx3g -jar $GATKHOME/GenomeAnalysisTK.jar -T SelectVariants -R data/ref_files/Parus_major_1.04.chrLGE22.fa -V vcf_files/samtools.chrLGE22.raw.snps.indels.vcf.gz -o vcf_files/samtools.chrLGE22.raw.snps.vcf.gz -selectType SNP -restrictAllelesTo BIALLELIC 
    
How many SNPs in the VCF files?

Extract and compare the number of SNPs called by each caller using bedtools. Bedtools is a very useful tools for working
with genome interval data in bed, VCF or GFF format.

    module load apps/gcc/5.2/bedtools
    bedtools intersect -a 
    
    
## Filtering SNPs

### Depth filters

### Base qualities, Strand Biases

### Region Filters

### GATK recommended Hard Filters

## ANGSD (SNP calling and Popgen analysis in low coverage data)

## SFS

Comparing sfs

