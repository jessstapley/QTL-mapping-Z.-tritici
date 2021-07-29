# RNAseq, DEG and topGO analysis 
This document contains a step-by-step account of how we analyzed the RNA seq data, performed a differential gene expression analysis and topGO analysis.  Here we have provided examples of code for each step. The code chunks will not run automatically "as is" it will need to be edited by the user.

List of the software used (not including all dependencies)
star v2.5.3a, cufflinks v2.1.1, R v3.6.0


## Check RNA sequence read quality with ```fastqc```
Example code

```
while read p; do
fastqc path_to_seq_data/RNAseq_${name}.fastq.gz -o fastqc_out/
done<RNAseq_sample.list

```
## Triming
We trimmed the raw RNA seq reads using ```trimmmomatic v0.35``
```
trimmomatic PE -threads 2 -phred33 ${name}_R1.fastq.gz ${name}_R2.fastq.gz ${out}/${name}.R1.trim.fq.gz ${out}/logs/${name}.R1.un.fq.gz ${out}/${name}.R2.trim.fq.gz ${out}/logs/${name}.R2.un.fq.gz ILLUMINACLIP:${adaptor}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 > ${out}/logs/${name}.trimmo.log  2> ${out}/logs/${name}.trimmo.err

```
## Generate genome files for ```star``` 
We ran STAR using the gff from Treindl et al 2021, accoding to the manual recommendations. As we useed a gff not a gtf we added the following option --sjdbGTFtagExonParentTranscript

```
fasta="/path_to_fasta/Ety1756_Epichloe_typhina_1756_33930528_v4.fna"
gtf="/path_to_gff/Epichloe_typhina.gff3"
GDIR="path_genome_files/ggIndex_Et"

STAR --runMode genomeGenerate --runThreadN 20 --genomeDir $GDIR  \
--genomeFastaFiles $fasta \
--sjdbGTFfile      $gtf \
--sjdbGTFtagExonParentTranscript Parent  \
--sjdbOverhang 124

```

## Map RNAseq reads to the reference genome using ```star```
```
GDIR="path_genome_files/ggIndex_Et"

STAR --runMode alignReads --runThreadN 20 --genomeDir $GDIR \
--readFilesIn /path_to_fastq/${name}_R1.fastq.gz /path_to_fastq/${name}_R2.fastq.gz \
--readFilesCommand zcat    \
--outFileNamePrefix ${name}_ \
--outSAMtype BAM SortedByCoordinate    \
--outSAMstrandField intronMotif 

```

## Extract the read counts per gene using R

We used R to count the number of reads per gene. This package requires a gtf file, so we used  ``cufflinks`` to convert the gff file to a gtf file.
```
gffread /path_to_gff/Epichloe_clarkii.gff3 -T -F -o /path_to_gtf/Epichloe_clarkii.gtf

```
Using R package  ``Subread`` to count the reads per feature (gene) and create an ouptut file called "Allcounts*"

```
names = read.table("bam.list")

file_names = as.vector(names[,1])
name.prefix <- gsub("_Aligned.sortedByCoord.out.bam", "", file_names)

fx = featureCounts(files = file_names, annot.ext="/oath_to_gft/Epichloe_typhina.gtf", isGTFAnnotationFile = TRUE, isPairedEnd = TRUE)

non_zero = rowSums(fx$counts)>0
counts = fx$counts[non_zero,]
colnames(counts) = paste0(name.prefix)

write.table(counts, file="Allcounts_Et", quote=FALSE)

```

## Identify differentially expressed genes using R

We used R package EdgeR to identify genes that were differentially expressed. We first plotted the data to check for outliers, using a clustering plot and a PCA. We used log counts to reduce the influence of genes with high counts and selected the top 500 most variable genes to speed up the process.

```
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("edgeR", "GO.db"))

library(edgeR)
library(gplots) 

counts.et = read.table("data/Allcounts_Et", header=T)
y <- DGEList(counts=counts.et, group=id.info.ord.et$treatment)

logcounts <- cpm(y,log=TRUE)
var_genes <- apply(logcounts, 1, var)

select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]

highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
```

We added treatment information to the variable counts data frame and then created two plots - a hierarchical clustering plot and a multidimensional scaling (MDS) plot. These plots ar useful to look for outliers or errors in your data.
```
colnames(highly_variable_lcpm) = c(paste0("plant",1:6), paste0("plate", 1:3))

pdf("plots/clustering_Et.pdf")
heatmap.2(highly_variable_lcpm, trace="none", main="Top 500 most variable genes across samples")
dev.off()

pdf("plots/MDS_Et.pdf")
plotMDS(y, col=c(rep("red",3), rep("blue", 4)))
dev.off()

```

## Differntial gene expression analysis

First we select the highly expressed genes, then created a DGEList . Variation in sequencing depth can influence variation in read counts, so we normalise for library size. Often the number of reads is used as the library size, but in this experiment the number of mapped reads is more appropiate because the samples from 'in planta' contain moslty plant RNA and the number of reads dos not reflect the sequencing depth at a gene.

```
library(edgeR)
library(GO.db) 
library(corrplot)
library(gplots) 

keep <- rowSums(cpm(counts.et)>10) >=2 
counts.et = counts.et[keep,]

et.dat = DGEList(counts=counts.et, genes=row.names(counts.et), lib.size = id.info.ord.et.rm$mapped, group = id.info.ord.et$treatment)

et.dat <- calcNormFactors(et.dat)
design <- model.matrix(~id.info.ord.et$treatment)
```
Then we estimate dispresion across all genes (Tags) and we perform a genewise exact test for differences in the means betwen groups. 

```
et.dat <- estimateDisp(et.dat)
et.dat <- estimateCommonDisp(et.dat)
et.dat <- estimateTagwiseDisp(et.dat)

# testing for differenitally expressed genes
xt.et.dat = exactTest(et.dat, pair=levels(et.dat$samples$group))
topTags(xt.et.dat) # most differentiated 

deg.et.dat <- decideTestsDGE(xt.et.dat, adjust.method="BH", p.value=0.05)
summary(deg.et.dat)

```
## GO analysis on DEGs

