# SNP calling and Filtering

We will use some whole genome data from 10 great tit individuals to look at SNP calling and filtering.
These 10 individuals were sampled in Europe and are a subset of the 29 birds that were sequenced and 
analysed in [Laine et al. (2016)](http://www.nature.com/ncomms/2016/160125/ncomms10474/full/ncomms10474.html).

We will focus on a small subset of the genome, calling SNPs on chrLGE22

The 10 birds were sequenced with 100bp paired-end Illumina reads to a mean depth of between 9-14x.
In the Laine et al. (2016) paper the SNP calling was performed using GATK, Platypus and ANGSD.

## Programs required
All the programs we will use are available on iceberg through the module system or from the genomics repository.

The following programs will be used:
    
* samtools v1.2 (available through the module system)
* bcftools v1.3 (available in the genomics repository)
* GATK v3.4 (available through the  module system)
* bedtools v2.5 (available through the module system)

## Obtaining the tutorial material

Type the following at the command terminal.

    git clone https://github.com/padraicc/uspopgen
    
This will download a folder called uspopgen into your current directory. Change directory into the uspopgen folder.
    
    cd uspopgen
    
Now take a look at the contents.

    ls
    
You will find folders called ```data/``` and ```scripts/``` and a ```README.md``` file containing the the text for this webpage.

## Data files

The data files we will use are located in the data folder

    ls data

**gvcf_files/** contains the 10 gVCF files needed for the GATK SNP calling.  
**bam_files/** cotains the BAM files needed for the samtools SNP calling.  
**ref_files/** contains the great tit reference genome files for chrLGE22.  

## SNP callers

We will call SNPs using two popular SNP calling programs: GATK and samtools/bcftools.

### GATK SNP calling

We will perform SNP calling of the 10 birds using the shell script called ```gatk_snp_calling.sh```. Due to time constraints
we will only run the last step of the GATK pipeline. Here we use the *GenotypeGVCFs* tool to produce the VCF files from the
10 gVCF files in ```data/gvcf_files```. (The gVCF used as input files were produced using the GATK *HaplotypeCaller* tool.) 

If you wish to take a look at the gatk command inside the shell script you can use ```less``` command as follows

    less -S scripts/gatk_snp_calling.sh

To run the script in the terminal just type the following:
    
    bash scripts/gatk_snp_calling.sh 
    
This will result in a compressed VCF file being written in the new *vcf_files* folder found in the current working 
directory. 

    ls vcf_files

### SAMtools and bcftools SNP calling

We will also use the samtools and bcftools programs to call SNPs from the BAM file for our 10 birds. The BAM contain the
alignments of the reads mapped to the great tit chrLGE22 reference genome. The bam files are located in the ```data/bam_files/```
folder. 

To run the samtools snp calling just type the following:

    bash scripts/samtools_snp_calling.sh
    
### Understanding the VCF format
Now that we have run the SNP calling programs we will take a look at the VCF format files contained in the ```vcf_files/```
folder.
 
The VCF format is composed of a header section where each line begins with '##' and the headers describing the columns
are located on the line starting with '#CHROM'. 

**Table 1** Description of VCF fields

| Column Number| Title | What it contains |  
|:--|:--|:--|  
| 1 | CHROM | Chromosome/scaffold of the reference genome |
| 2 | POS | Position on the chromosome/scaffold given in CHROM (1-based) |
| 3 | ID | ID for the SNP/INDEL (Here always '.')  |
| 4 | REF | Allele in the reference genome  |
| 5 | ALT | The variant allele |
| 6 | QUAL | The variant quality score (phred-scale) |
| 7 | FILTER | The field storing filtering information as a string ('.' in an unfiltered file) |
| 8 | INFO | List of variant annotations delimited by ';'  |
| 9 | FORMAT | Describes the format for the sample genotype information in column 10 onward |
| 10 to end | Sample id (from SM tag in BAM file RG string) | Sample genotype information (Usually also stores the genotype qualities (GQ) and genotype likelihoods (GL or PL (phred-scaled) |
 
 
The INFO field contains a lot of annotations for the variant site (e.g Depth, Mapping Qualtity etc.) and
these values may be used for filtering. The format for the genotype information is explained in the header
of VCF file. Each row following the header section is a variant site in the VCF, either a SNP or an INDEL. For more
information on the VCF format see [here.](http://samtools.github.io/hts-specs/VCFv4.2.pdf) and [here.](http://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it)

Now, take a look at the GATK and samtools vcf files.

    less -S vcf_files/gatk.chrLGE22.raw.snps.indels.vcf.gz
    less -S vcf_files/samtools.chrLGE22.raw.snps.indels.vcf.gz 

If you want to exclude the header section, try the following.
    
    zgrep -v ^## vcf_files/gatk.chrLGE22.raw.snps.indels.vcf.gz | less -S

To find an explanation of the info field in the samtools VCF, look at the last 25 lines of the header section
    
    zgrep ^# vcf_files/samtools.chrLGE22.raw.snps.indels.vcf.gz | tail -n 25
    
Note that different SNP callers will produce some different annotations in the INFO field and genotype
fields. 

**Q1.** What differences do you see between the INFO field of the samtools VCF and the GATK VCF files?

Samtools is less well documented with regard to exactly what some of its annotations are, beyond what can be understood from the VCF header. 
GATK has a much richer documentations which can be browsed [here.](https://www.broadinstitute.org/gatk/guide/tooldocs/#VariantAnnotations)
Click on Variant annotations for more info.


### Comparing the output from the two callers

We will count and compare the unfiltered SNP calls using some common command line tools and the bedtools program.

Count the number of variants in the vcf file using zgrep (grep for compressed files).

    zgrep -cv ^# vcf_files/gatk.chrLGE22.raw.snps.indels.vcf.gz 
    zgrep -cv ^# vcf_files/samtools.chrLGE22.raw.snps.indels.vcf.gz 
    
**Q2.** Which tool calls the most variants?

Now we will extract only biallelic SNPs and discard the INDELs from the VCF using GATK's SelectVariants tool.
    
    module load apps/binapps/GATK
    java  -Xmx3g -jar $GATKHOME/GenomeAnalysisTK.jar -T SelectVariants -R data/ref_files/Parus_major_1.04.chrLGE22.fa -V vcf_files/gatk.chrLGE22.raw.snps.indels.vcf.gz -o vcf_files/gatk.chrLGE22.raw.snps.vcf.gz -selectType SNP -restrictAllelesTo BIALLELIC 

We will also extract the SNPs from the samtools vcf, but first we need to index it, as the GATK tools require this.

    tabix -p vcf vcf_files/samtools.chrLGE22.raw.snps.indels.vcf.gz
    java  -Xmx3g -jar $GATKHOME/GenomeAnalysisTK.jar -T SelectVariants -R data/ref_files/Parus_major_1.04.chrLGE22.fa -V vcf_files/samtools.chrLGE22.raw.snps.indels.vcf.gz -o vcf_files/samtools.chrLGE22.raw.snps.vcf.gz -selectType SNP -restrictAllelesTo BIALLELIC 
    
**Q3.** How many SNPs in the VCF files?

Extract and compare the number of SNPs called by each caller using bedtools. Bedtools is a very useful tools for working
with genome interval data in bed, VCF or GFF format.

    module load apps/gcc/5.2/bedtools
    
We will first extract the SNPs that were called by both callers.
    
    bedtools intersect -header -a vcf_files/gatk.chrLGE22.raw.snps.vcf.gz -b vcf_files/samtools.chrLGE22.raw.snps.vcf.gz | bgzip > vcf_files/both.chrLGE22.raw.snps.vcf.gz
    
Next we will find SNPs called by only one of the callers.    
    
    bedtools subtract -header -a vcf_files/gatk.chrLGE22.raw.snps.vcf.gz -b vcf_files/samtools.chrLGE22.raw.snps.vcf.gz | bgzip > vcf_files/gatk_only.chrLGE22.raw.snps.vcf.gz
    bedtools subtract -header -a vcf_files/samtools.chrLGE22.raw.snps.vcf.gz -b vcf_files/gatk.chrLGE22.raw.snps.vcf.gz | bgzip >  vcf_files/samtools_only.chrLGE22.raw.snps.vcf.gz
    
**Q4.** How many SNPs did both callers call? How many by only one of the callers?

## Filtering SNPs

We now have the raw SNPs in VCF format. However, these files are likely to contain a number of false positives, so we need
to apply some hard filters to remove sites that may not be true SNPs. Hard filtering refers to setting a threshold on
a number of SNP properties (many of those in the listed INFO field of the VCF) and removing any SNPs that fail to meet 
that threshold.

Some common filters applied to SNP calls include:

* Depth filters (Setting a minimum and maximum depth for a site)  
* Minimum SNP quality (Requiring a minimum quality in the QUAL field of the SNP)  
* Minimum RMS mapping quality for SNPs    
* Strand Bias (Filtering sites where the number of reference and non-reference reads are highly correlated with the strands of the reads)  
* Hardy-Weinberg Equilibrium (Filter sites that show significant deviation from Hardy-Weinberg Expectations) 

For a summary of SNP filtering applied to whole genome resequencing studies in birds see [here.](https://www.dropbox.com/s/xa0bndtz42i1uft/snp_filtering_avian_studies.pdf?dl=0)

One of the simplest thresholds that can be applied is a minimum quality score. Here we will use the a *bcftools filter*
to filter our samtools SNP VCF. Documentation on this tool can be found [here.](https://samtools.github.io/bcftools/bcftools.html#filter)

    export PATH=/usr/local/extras/Genomics/apps/bcftools/1.3/bin/:$PATH
    bcftools filter -e "QUAL<30" -s "LOW_QUAL" vcf_files/samtools.chrLGE22.raw.snps.vcf.gz -O z -o vcf_files/samtools.chrLGE22.qual_filtered.snps.vcf.gz

**Q5.** How many SNPs PASS at the QUAL < 30 filter for samtools VCF?

(**Note** that the number of SNPs in the QUAL filtered samtools VCF is now more similar to the GATK raw VCF. This is because the 
GenotypeGVCF tool by default excludes calls with QUAL less than 30. See the documentation on GenotypeGVCF tools 
[here](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php))

We can also apply a depth filters with ```bcftools filter```. (Note that the '||' stand for the logical 'or'). The mean depth
for sites was 106 across samples. Here we apply a filter on sites with greater than twice the mean or less than half the mean depth.
    
    bcftools filter -e "QUAL<30 || MAX(DP)>213 || MIN(DP)<53 " -s "LOW_QUAL" vcf_files/samtools.chrLGE22.raw.snps.vcf.gz -O z -o vcf_files/samtools.chrLGE22.qual_dp_filtered.snps.vcf.gz
 
Count the number of filtered sites as follows.
 
    zgrep -v ^# vcf_files/samtools.chrLGE22.qual_dp_filtered.snps.vcf.gz | grep -cw LOW_QUAL
    
 **Q6.** How many SNPs are excluded when we apply both the depth and quality filter to the samtools VCF?   
     
### GATK recommended Hard Filters

In GATK there are two options for filtering SNPs, either hard or soft filtering. The first is to use a 'soft' filtering 
approach by using known variants to carry out Variant Quality Score Recalibration. This approach requires a large number 
(100,000s) of known SNPs and may not be useful for researchers working on non-model organisms, or species without 
a large existing set of polymorphism data. Also, this approach may not be used for targeted resequencing such 
GBS or RADseq data.

In cases where VQSR can not be used, the GATK developers recommend a set of hard filters for filtering SNPs. Further details
on these recommended filters can be seen [here](http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set)

To apply the hard filters run the following command on the GATK VCF containing only SNPs.

    java -Xmx3g -jar $GATKHOME/GenomeAnalysisTK.jar -T VariantFiltration -R data/ref_files/Parus_major_1.04.chrLGE22.fa -V vcf_files/gatk.chrLGE22.raw.snps.vcf.gz --filterExpression "QD<2.0||FS>60.0||MQ<40.0||MQRankSum<-12.5||ReadPosRankSum<-8.0" --filterName "GATK_hard_snp_filter" -o vcf_files/gatk.chrLGE22.hard_filtered.snps.vcf.gz

(**Note** that the VariantFiltration tool issues a warning at sites the MQRankSum and ReadPosRankSum annotations.
See why this is the case [here](http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set))

Count the number of SNPs that PASS the filters.

    zgrep -v ^# vcf_files/gatk.chrLGE22.hard_filtered.snps.vcf.gz | grep -cw PASS

**Q7.** How many SNPs were filtered using the GATK recommended filters? Do the recommended GATK hardfilters filter based on th QUAL score?


The GATK recommendation may need tweaking. It is worthwhile to extract the annotations from the VCF file and plot them to get
an idea where to set the cut offs for the various annotations. A tool such as the VariantsToTable (see [here](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php)) 
tools is a convenient way to extract the annotations to a text file that may be read into R.


## Region Filters

Finally, if you have time you can look at how to filter SNPs that fall in repetitive regions of the genome.

Another common additional filter applied to filtering SNPs in whole genome data is to exclude SNPs that fall in 
repetitive regions of the genomes where mapping may be problematic. Here we will further filter the GATK hard filtered 
VCF file to exclude SNPs in repetitive regions using a bed file located at ```data/ref_files/chrLGE22.reps.bed``` which 
specifies the repetitive regions of chrLGE22.

**Q8.** How might you exclude SNPs in our filtered GATK VCF file within repetitive regions using bedtools?


*Solutions to all the questions are available from here* 