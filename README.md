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
in the format fields. 

Q1. What differences do you see between the INFO field of the samtools VCF and the GATK VCF files?

### Comparing the output from the two callers

We will count and compare the unfiltered SNP calls using some common command line tools and the bedtools program.

Count the number of variants in the vcf file using zgrep (grep for compressed files).

    zgrep -cv ^# vcf_files/gatk.chrLGE22.raw.snps.indels.vcf.gz 
    zgrep -cv ^# vcf_files/samtools.chrLGE22.raw.snps.indels.vcf.gz 
    
Q2. Which tool calls the most variants?

Now we will extract only biallelic SNPs and discard the INDELs from the VCF using GATK's SelectVariants tool
    
    module load apps/binapps/GATK
    java  -Xmx3g -jar $GATKHOME/GenomeAnalysisTK.jar -T SelectVariants -R data/ref_files/Parus_major_1.04.chrLGE22.fa -V vcf_files/gatk.chrLGE22.raw.snps.indels.vcf.gz -o vcf_files/gatk.chrLGE22.raw.snps.vcf.gz -selectType SNP -restrictAllelesTo BIALLELIC 

We will extract the SNPs from the samtools vcf, but first we need to index it, as the GATK tools require this.

    tabix -p vcf vcf_files/samtools.chrLGE22.raw.snps.indels.vcf.gz
    java  -Xmx3g -jar $GATKHOME/GenomeAnalysisTK.jar -T SelectVariants -R data/ref_files/Parus_major_1.04.chrLGE22.fa -V vcf_files/samtools.chrLGE22.raw.snps.indels.vcf.gz -o vcf_files/samtools.chrLGE22.raw.snps.vcf.gz -selectType SNP -restrictAllelesTo BIALLELIC 
    
Q3. How many SNPs in the VCF files?

Extract and compare the number of SNPs called by each caller using bedtools. Bedtools is a very useful tools for working
with genome interval data in bed, VCF or GFF format.

    module load apps/gcc/5.2/bedtools
    bedtools intersect -header -a vcf_files/gatk.chrLGE22.raw.snps.vcf.gz -b vcf_files/samtools.chrLGE22.raw.snps.vcf.gz | bgzip > vcf_files/both.chrLGE22.raw.snps.vcf.gz
    
Now we will find SNPs called by only one of the callers.    
    
    bedtools subtract -header -a vcf_files/gatk.chrLGE22.raw.snps.vcf.gz -b vcf_files/samtools.chrLGE22.raw.snps.vcf.gz | bgzip > vcf_files/both.chrLGE22.raw.snps.vcf.gz
    bedtools subtract -header -a vcf_files/samtools.chrLGE22.raw.snps.vcf.gz -b vcff_files/gatk.chrLGE22.raw.snps.vcf.gz | bgzip >  vcf_files/both.chrLGE22.raw.snps.vcf.gz
    
Q4. How many SNPs were called by both callers? How many by only one of the callers?

## Filtering SNPs

We now have the raw SNPs in VCF format. However, these files are likely to contain a number of false positives, so we need
to apply some hard filters to remove sites that may not be true SNPs. Hard filtering refers to setting a threshold on
a number of SNP properties (many of those in the listed INFO field of the VCF) and removing any SNPs that fail to meet a
particular meet that threshold.

A number of common filters applied to SNP calls include
    * Depth filters (setting a minimum and maximum depth for a site)
    * Minimum SNP quality (Requiring a minimum quality in the QUAL field of the SNP)
    * Minimum RMS mapping quality for SNPs
    * Allele Balance (filtering sites where the fraction of non-reference reads is too low) 
    * Strand Bias 

For a summary of SNP filtering applied to whole genome resequencing studies in birds see [here](https://www.dropbox.com/s/xa0bndtz42i1uft/snp_filtering_avian_studies.pdf?dl=0)

One of the simplest thresholds that can be applied is a minimum quality score. Here we will apply this to samtools VCF
containing only SNPs. Here we will 



### GATK recommended Hard Filters

In GATK there are two options for filtering SNPs, either hard of soft filtering. The first is to use a 'soft' filtering 
approach by using known variants to carry out Variant Quality Score Recalibration. This approach requires a large number 
(100,000s) of known SNPs and maynot be useful for researcher working on non-model organisms, or species without 
large existing set of polymorphism data. Also, this approach may not be used for targeted reseuqencing such 
GBS or RADseq data.

In cases where VQSR can not be used, the GATK developer recommend a set of hard filters for filtering SNPs. Further details
on these recommended filters can be seen here [here](http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set)

To apply the hard filters run the following command on the GATK VCF containing only SNPs.

    java -jar GenomeAnalysisTK.jar -T VariantFiltration -R data/ref_files/Parus_major_1.04.chrLGE22.fa -V vcf_files/gatk.chrLGE22.raw.snps.vcf.gz --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"  --filterName "GATK_hard_snp_filter" -o vcf_files/gatk.chrLGE22.hard_filtered.snps.vcf.gz

Count the number of SNPs that PASS the filters

    zgrep -v ^# PASS vcf_files/gatk.chrLGE22.hard_filtered.snps.vcf.gz

## Region Filters

Another common additional filter applied to filtering SNPs in whole genome data is to exclude SNPs that fall in 
repetitve regions of the genomes. Here we will further filter the GATK hard filtered VCF file to exclude SNPs in repetive
regions using a bed file which specifies the repetitive regions of chrLGE22

How might you exclude SNPs in our GATK VCF file within repetitive regions using bedtools subtract?


*Solutions to all the questions can be downloaded from here* 