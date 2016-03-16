# SNP calling and Filtering

We will use some whole genome data from 10 great tit idnviduals to look at SNP calling and filtering.
These 10 individuals were sampled in Europe and area subset of the 29 birds that were sequenced and 
analysed in [Laine et al. (2016)](http://www.nature.com/ncomms/2016/160125/ncomms10474/full/ncomms10474.html).

We will foucs on a small subset of the genome, calling SNPs on chrLGE22

The 10 birds were sequenced with 100bp paired-end Illumina reads to a mean depth of between 9-14x.
In the Laine et al. (2016) paper the SNP calling was performed using GATK, Platypus and ANGSD.

## Programs required
All the programs we will use are available on iceberg through the module system or from the genomics repository.

The following programs will be used:
    
* samtools v1.2 (available through the module system)
* bcftools v1.3 (available in the genomics repository)
* GATK v3.4 (available through the  module system)
* ANGSD v0.911

## Obtaining the tutorial material

Type the following at the command terminal

    git clone https://github.com/padraicc/uspopgen

Change directory into the uspopgen folder
    
    cd uspopgen
    ls
    
You will find folders called data/ and scripts/ and a README.md files containing the the text for this webpage.

## Data files

The data files we will use are located in the data folder

    ls data

**gvcf_files/** contains the 10 gVCF files needed for the GATK SNP calling.  
**bam_files/** cotains the BAM files for samtools SNP calling.  
**ref_files/** contains the reference genome file for chrLGE22.  

## SNP callers

We will call SNPs using GATK and samtools/bcftools.

### GATK SNP calling

We will perform SNP calling of the 10 birds using the shell script called **gatk_snp_calling.sh**. Due to time constraints
we will run the last strep of the GATK pipeline. Here we use the *GenotypeGVCFs* tool to produce the VCF files from the
10 gVCF files in */data/gvcf_files*. (The gVCF used as input files were produced using the GATK HaplotypeCaller program.) 

To run the script in the terminal just type the following:
    
    bash scripts/gatk_snp_calling.sh 

This will result in a compressed VCF file being written in the new *vcf_files* folder found in the current working 
directory.

    ls vcf_files

### SAMTOOLS/BCFTOOLS SNP calling

    bash scripts/samtools_snp_calling.sh
    
### Understanding the VCF format
First, we will take a look at the VCF format, as this is the output format for most SNP calling and genotyping programs
 
The VCF format is composed of a header section where each line begins with '##' and the headers describing the columns
are located on the line starting with '#chrom'. The columns of the VCF from left to right are the chromosome/scaffold,  
position within the reference chromosome/scaffold, an ID (here always ‘.’), the reference allele, the variant allele,
the SNP/INDEL quality score, the filter field, info field, the sample format field and the samples make up the remaining
columns (10 in this case).

The INFO field contains a lot of annotations for the variant site (e.g Depth, Mapping Qualtity etc.) and
these values may be used for filtering (see below). The format for the genotype information is explained in the header
of VCF file. Each row following the header section is a variant site in the VCF, either a SNP or an indel.

Now, take a look at the GATK and samtools vcf files.

    less -S vcf_files/gatk.chrLGE22.raw.snps.indels.vcf.gz
    less -S vcf_files/samtools.chrLGE22.raw.snps.indels.vcf.gz 

If you want to exclude the header section, try the following.
    
    zgrep -v ^## vcf_files/gatk.chrLGE22.raw.snps.indels.vcf.gz | less -S
    
Note that different SNP callers will have some differences in the annotations present in the INFO field and differences
in the format fields. What differences do you see between the INFO field of the samtools VCF and the GATK VCF files?

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
    bedtools intersect -header -a vcf_files/gatk.chrLGE22.raw.snps.vcf.gz -b vcf_files/samtools.chrLGE22.raw.snps.vcf.gz > vcf_files/both.chrLGE22.raw.snps.vcf.gz
    
    
## Filtering SNPs

### Depth filters

### Base qualities, Strand Biases

### Region Filters

### GATK recommended Hard Filters

## ANGSD (SNP calling and Popgen analysis in low coverage data)

## SFS

Comparing sfs

