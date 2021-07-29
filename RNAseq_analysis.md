# RNAseq, DEG and topGO analysis 
This document contains a step-by-step account of how I analyzed the RNA seq data and performed a differential gene expression analysis.  Here I have provided examples of code for each step. The code chunks will not run automatically "as is" it will need to be edited by the user.

## Check RNA sequence read quality with ```fastqc```
Example code

```
while read p; do
fastqc path_to_seq_data/RNAseq_${name}.fastq.gz -o fastqc_out/
done<RNAseq_sample.list

```
## Triming
I trimmed the raw RNA seq reads using ```trimmmomatic v0.35``
```
trimmomatic PE -threads 2 -phred33 ${name}_R1.fastq.gz ${name}_R2.fastq.gz ${out}/${name}.R1.trim.fq.gz ${out}/logs/${name}.R1.un.fq.gz ${out}/${name}.R2.trim.fq.gz ${out}/logs/${name}.R2.un.fq.gz ILLUMINACLIP:${adaptor}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 > ${out}/logs/${name}.trimmo.log  2> ${out}/logs/${name}.trimmo.err

```
## Generate genome files for ```star v2.5.3a``` 
I ran STAR using accoding to the manual recommendations. When using a gff not a gtf you need to add the following option --sjdbGTFtagExonParentTranscript

```
fasta="/path_to_fasta/reference.fa"
gtf="/path_to_gff/reference.gff3"
GDIR="path_genome_files/ggIndex"

STAR --runMode genomeGenerate --runThreadN 20 --genomeDir $GDIR  \
--genomeFastaFiles $fasta \
--sjdbGTFfile      $gtf \
--sjdbGTFtagExonParentTranscript Parent  \
--sjdbOverhang 124

```

## Map RNAseq reads to the reference genome using ```star```
```
GDIR="path_genome_files/ggIndex"

STAR --runMode alignReads --runThreadN 20 --genomeDir $GDIR \
--readFilesIn /path_to_fastq/${name}_R1.fastq.gz /path_to_fastq/${name}_R2.fastq.gz \
--readFilesCommand zcat    \
--outFileNamePrefix ${name}_ \
--outSAMtype BAM SortedByCoordinate    \
--outSAMstrandField intronMotif 

```

## Extract the read counts per gene using R

I used R (v3.6.0) package Rsubread to count the reads per feature (gene) and create an ouptut file called "Allcounts"

```
names = read.table("bam.list")

file_names = as.vector(names[,1])
name.prefix <- gsub("_Aligned.sortedByCoord.out.bam", "", file_names)

fx = featureCounts(files = file_names, annot.ext="/path_to_gft/reference.gtf", isGTFAnnotationFile = TRUE, isPairedEnd = TRUE)

non_zero = rowSums(fx$counts)>0
counts = fx$counts[non_zero,]
colnames(counts) = paste0(name.prefix)

write.table(counts, file="Allcounts", quote=FALSE)

```

## Identify differentially expressed genes using R

I used R (v3.6.0) package EdgeR to identify genes that were differentially expressed between strains. I first plotted the data to check for outliers, using a clustering plot and a PCA. 

```
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("edgeR", "GO.db"))

library(edgeR)
library(gplots) 

counts = read.table("data/Allcounts", header=T)
y <- DGEList(counts=counts, group=id.info$treatment)

logcounts <- cpm(y,log=TRUE)
var_genes <- apply(logcounts, 1, var)

select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]

highly_variable_lcpm <- logcounts[select_var,]

pdf("plots/clustering.pdf")
heatmap.2(highly_variable_lcpm, trace="none", main="Top 500 most variable genes across samples")
dev.off()

pdf("plots/MDS.pdf")
plotMDS(y, col=c(rep("red",3), rep("blue", 4)))
dev.off()

```

## Differntial gene expression analysis

First I select the highly expressed genes, then created a DGEList. Variation in sequencing depth can influence variation in read counts, so we normalise for library size (number of reads). For the invivio samples I used the number of mapped reads instead of the total number of reads becasue many reads came from plant DNA and did not map to the Z. tritici genome, so total number of reads did not correspond to sequencing depth. 

```
library(edgeR)

keep <- rowSums(cpm(counts)>10) >=2 
counts = counts[keep,]

dat = DGEList(counts=counts, genes=row.names(counts), lib.size = id.info$mapped, group = id.info$treatment)

dat <- calcNormFactors(dat)
design <- model.matrix(~id.info$treatment)
```
Then I estimated dispresion across all genes (Tags) and performed a genewise exact test for differences in the means betwen groups (strains). 

```
dat <- estimateDisp(dat)
dat <- estimateCommonDisp(dat)
dat <- estimateTagwiseDisp(dat)

# testing for differenitally expressed genes
xt.dat = exactTest(dat, pair=levels(dat$samples$group))
topTags(xt.dat) # most differentiated 

deg.dat <- decideTestsDGE(xt.dat, adjust.method="BH", p.value=0.05)
summary(deg.dat)

```

