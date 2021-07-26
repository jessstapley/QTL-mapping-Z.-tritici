# Variant calling with GATK Short Variant Discovery

## Finding variants in an individual sample
We first split each genome into 100000bp regions. We then run ```HaplotypeCaller``` on each region, looping thorugh all regions to produce g.vcf file for each region. We then join the regions with picard tool ```GatherVcfs```.

```
for i in {1..nRegions}
do 
java -jar /path_to_package/gatk-package-4.1.2.0-local.jar HaplotypeCaller \
-R /path_to_ref_genome/Ref.fa \
--sample-ploidy 1 \
-I ${data}/${name}_sort_10_dup.bam \
-L /path_to_regions_file/regions_100000_${i}.intervals \
-O ${path_out}/${name}/${name}_r${i}.raw.snps.indels.g.vcf

picard SortVcf INPUT=${path_out}/${name}/${name}_r${i}.raw.snps.indels.g.vcf \
OUTPUT=${path_out}/${name}/${name}_r${i}_sort.raw.snps.indels.g.vcf 
done	

for i in {1..nRegions}
do 
picard GatherVcfs INPUT=${path_out}/${name}/${name}_r${i}_sort.raw.snps.indels.g.vcf \
OUTPUT=${path_out}/${name}/${name}_sort.raw.snps.indels.g.vcf 
done

```

These read Filters are automatically applied to the data by GATK before processing by HaplotypeCaller.

HCMappingQualityFilter, minQual =20
MalformedReadFilter - avoids crashing on malformed reads
BadCigarFilter
UnmappedReadFilter - only mapped reads
NotPrimaryAlignmentFilter -This filter recognizes the SAM flag that identifies secondary alignments. It is intended to ensure that only records that are likely to be mapped in the right place, and therefore to be informative, will be used in analysis. To be clear, it does NOT filter out read records that are supplementary alignments.
FailsVendorQualityCheckFilter - This filter recognizes the SAM flag corresponding to the vendor quality check.
DuplicateReadFilter - This filter recognizes the SAM flag set by MarkDuplicates.
MappingQualityUnavailableFilter - This filter is intended to ensure that only reads that are likely to be mapped in the right place, and therefore to be informative, will be used in analysis.


## Combine individual gvcf files
The GVCF files for each invidivual from each cross was combined using ``` CombineGVCFs```

```
java -jar /path_to_package/gatk-package-4.1.2.0-local.jar CombineGVCFs \
-R /path_to_reference/Ref.fna \
--variant /path_to_data/GATK_HC/${name}/${name}_sort.raw.snps.indels.g.vcf \
...n Samples
-O /result_dir/Cross1A5_Combined_sort.raw.snps.indels.g.vcf

```

## Joint Calling to find SNPs
Next  we used the GATK joint genoytoing function to identify SNPs across all samples.

```
java -jar /path_to_package/gatk-package-4.1.2.0-local.jar GenotypeGVCFs \
-R /path_to_reference/Ref.fna \
--sample-ploidy 1 \
--variant /path_to_data/Cross1A5_Combined_sort.raw.snps.indels.g.vcf \
-O /result_dir/Cross1A5_JC_sort.raw.snps.indels.g.vcf

```

