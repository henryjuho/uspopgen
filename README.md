# SNP calling and Filtering

We will use some whole genome data from the great tit to look at SNP calling and filtering.
The data we will 10 great tit individuals from Europe, they are a subset of the 29 birds that were sequenced and 
in[Laine et al. (2016)](http://www.nature.com/ncomms/2016/160125/ncomms10474/full/ncomms10474.html).

The 10 birds were sequenced with illumina paired-end data with 100bp reads to a mean depth of between 9-14x.
In the Laine et al. (2016) paper the SNP calling was performed using GATK, Platypus and ANGSD.

## Programs required
All the programs we will use are available on iceberg through the module system or from the genomics repository.

The following programs will be used:
    
* samtools v1.2 (available through the module system)
* bcftools v1.3 (available in the genomics repository)
* GATK v3.4 (available through the  module system)
* ANGSD v

## SNP callers

We will call SNPs using GATK and samtools/bcftools.

### Understanding the VCF format
First, we will take a look at the VCF format, as this is the output format for most SNP calling and genotyping programs
 

### SAMTOOL/BCFTOOLS

    bash scripts/samtools_snp_calling.sh


### GATK

    bash scripts/gatk_snp_calling.sh 
   
### Comparing the output from different callers

We will count and compare the unfiltered SNP calls using some common command line tools and the bedtools program.

Count the number of variants in the vcf file using zgrep (grep for compressed files).

    zgrep -cv ^# vcf_files/gatk.chrLGE22.raw.snps.indels.vcf.gz 
    zgrep -cv ^# vcf_files/samtools.chrLGE22.raw.snps.indels.vcf.gz 
    
Which tool calls the most variants?

Now we will extract only SNPs and discard the INDELs from the VCF using GATK's SelectVariants toll
    
    module load apps/binapps/GATK
    
    
    
Count the number of SNPs in the VCF files.

Extract and compare the number of SNPs called by each caller using bedtools. Bedtools is a very usefule tools for working
with genome interval data in bed, VCF or GFF format.

    module load apps/gcc/5.2/bedtools
    
    
More information on these bedtools may be found here




## Filtering SNPs

Link to file on filtering used in population genomic studies in birds

### Depth filters

### Base qualities, Strand Biases

### Region Filters

### GATK recommended Hard Filters



## ANGSD (SNP calling and Popgen analysis in low coverage data)

## SFS

Comparing sfs

