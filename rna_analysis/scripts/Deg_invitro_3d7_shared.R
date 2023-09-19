# DEG mapped to 3D7

# need to replace PATH with path to working directory !!

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
BiocManager::install("GO.db")

library(edgeR)
library(DESeq2)
library(GO.db) 

library(corrplot)
library(gplots) 

counts = read.table("PATH/data/counts_invitro_3d7.txt", header=T)
head(counts)

id.counts = as.data.frame(names(counts))

id.info <- read.table("PATH/data/sample_info_invitro.txt", header=TRUE)
head(id.info)
names(id.info)[1] <- "id"
id.info$id # sample IDs - 

# add library size - number of reads
g3d7 <- read.table("PATH/data/mapStats_3d7.txt")
head(g3d7)
names(g3d7) <- c("id", "total", "mapped")
g3d7$pc.mapped <- g3d7$mapped/g3d7$total

head(g3d7)
head(id.info)
id.info <- merge(id.info, g3d7)

names(id.counts)[1] = "id"
dim(counts) # 11701 27
dim(id.info) # 27
id.infor.ordered = merge(id.counts, id.info, all.x=TRUE, sort=FALSE)
head(id.infor.ordered)
head(counts)

pdf("raw_counts.pdf")
# plot counts of samples belonging to the same group
plot(counts[,4], counts[,6], main = "raw counts", xlab = "B2", ylab = "B3")
abline(0,1)
# plot counts of samples belonging to different groups
plot(counts[,5], counts[,1], main = "raw counts", xlab = "B2", ylab = "A1")
abline(0,1)

# correlation between pairs of samples
corrplot(cor(counts), method="color")
dev.off()
# the correlation plot should show a higher correlation between the samples of the same group compared to that one between samples of diffrent groups.
#### you should normalise first becasue the genes with a lot of expression will flood the data


# we create a DGEList element
head(id.infor.ordered)
# group could be strain or treatment.. first remove mutants
# 3 SRR7076853 
# 7 SRR7076858 
# 27 SRR7076861

names(counts)[3]
names(counts)[7]
names(counts)[27]
dim(counts)# 11701    27

counts.all <- counts
counts <- counts.all[, -c(3, 7, 27)]
dim(counts)# 11701    24

#######
head(id.infor.ordered)
id.infor.ordered.all <- id.infor.ordered
id.infor.ordered <- id.infor.ordered.all[-c(3, 7, 27),]
dim(id.infor.ordered)

id.infor.ordered[,c(1,4)]
names(counts) <- id.infor.ordered$LibraryName


# split 3D versus 1A
id.infor.ordered$order <- seq(1, length(id.infor.ordered$id))
id.infor.3d <- subset(id.infor.ordered, id.infor.ordered$STRAIN=="ST99CH_3D1" | id.infor.ordered$STRAIN=="ST99CH_3D7")
counts.3d <- counts[,c(id.infor.3d$order)]
dim(counts.3d) # 12

write.csv(counts.3d, file="counts_invitro_g3d_3d13d7.csv")

id.infor.1a <- subset(id.infor.ordered, id.infor.ordered$STRAIN=="ST99CH_1A5" | id.infor.ordered$STRAIN=="ST99CH_1E4")
counts.1a <- counts[,c(id.infor.1a$order)]
dim(counts.1a) # 12
head(counts.1a)
write.csv(counts.1a, file="counts_invitro_g3d_1a5x1e4.csv")

# split into treatments and compare between strains
head(counts.3d)
id.infor.3d$order <- seq(1, length(id.infor.3d$id))

id.infor.3d.mm <- subset(id.infor.3d, id.infor.3d$trt=="MM")
counts.3d.mm <- counts.3d[,c(id.infor.3d.mm$order)]

id.infor.3d.ysb <- subset(id.infor.3d, id.infor.3d$trt=="YSB")
counts.3d.ysb <- counts.3d[,c(id.infor.3d.ysb$order)]
head(counts.3d.ysb)
head(counts.3d.mm)

####### mm 3d
y <- DGEList(counts=counts.3d.mm, group=id.infor.3d.mm$STRAIN)
y

logcounts <- cpm(y,log=TRUE)
# we use log counts to diminish the importance of high number counts.

dim(logcounts)
# to speed up the process, we select the 500 most variable genes.

var_genes <- apply(logcounts, 1, var)
head(var_genes)

# we select the 500 most variable log-cpm
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)

# Subset logcounts matrix (only keep the 500 most variables log-cpm)
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

colnames(highly_variable_lcpm) = id.infor.3d.mm$STRAIN

###### exploratory plots #############

pdf("clustering.pdf")
# we apply a hierarchical clustering on the samples and on the rows.
heatmap.2(highly_variable_lcpm, trace="none", main="Top 500 most variable genes across samples")
dev.off()
# splitting by strainn

pdf("MDS.pdf")
# Multidimensional scaling (MDS) plot:
plotMDS(y, col=c("red", "blue", "red", "red", "blue", "blue"))
# remember that MDS (and PCA) plots are useful yet partial representations of the data.
dev.off()

# alot of variation in 3D7, virtually none in 3d1



################ Differential expression analysis
# using highly expressed genes and comparing between heat and control


#Colcounts.CH.group - this file has id, treatment and reads information, you will need something similar. 
# how many reads (library size) did you get for each individual   

# select only highly expressed genes 
head(counts.3d.mm)
keep <- rowSums(cpm(counts.3d.mm)>10) >=2 
Hcounts.3d.mm = counts.3d.mm[keep,]
dim(Hcounts.3d.mm) #9492

# create a DGEList element

mm3d.str = DGEList(counts=Hcounts.3d.mm, genes=row.names(Hcounts.3d.mm), lib.size = id.infor.3d.mm$total, group = id.infor.3d.mm$STRAIN)
str(mm3d.str)
# exploratory plot - Check for outliers. ideally individuals of the same group should group together, but dont worry too much if they done. Good 
plotMDS(mm3d.str, col=c("red", "blue", "red", "red", "blue", "blue"))

# calculate normalisation factors, normalising on library size
mm3d.str <- calcNormFactors(mm3d.str)
design <- model.matrix(~id.infor.3d.mm$STRAIN)

# estimate disperson
mm3d.str <- estimateDisp(mm3d.str)
mm3d.str <- estimateCommonDisp(mm3d.str)
mm3d.str <- estimateTagwiseDisp(mm3d.str)

# testing for differenitally expressed genes
et.mm3d.str = exactTest(mm3d.str, pair=levels(mm3d.str$samples$group))
topTags(et.mm3d.str) # most differentiated 

dim(et.mm3d.str) # 8775 3
deg.mm3d.str <- decideTestsDGE(et.mm3d.str, adjust.method="BH", p.value=0.05)
summary(deg.mm3d.str)
#           ST99CH_3D7-ST99CH_3D1
# Down                    3512
# NotSig                  2658
# Up                      3322

# more differentiated than not!! 
3512+3322
sigTags.mm3d.str = topTags(et.mm3d.str, n=6834)
# write out the list of DEGs

write.table(sigTags.mm3d.str, file="PATH/out/deg_mm3d7v3d1_high10sig.txt", quote=FALSE, row.names=FALSE)
write.table(et.mm3d.str, file="PATH/out/eT_mm3d7v3d1_high10.txt", quote=FALSE, row.names=FALSE)


####### YSB 3d
y <- DGEList(counts=counts.3d.ysb, group=id.infor.3d.ysb$STRAIN)
y
logcounts <- cpm(y,log=TRUE)
# we use log counts to diminish the importance of high number counts
# to speed up the process, we select the 500 most variable genes.

var_genes <- apply(logcounts, 1, var)
# we select the 500 most variable log-cpm
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)

# Subset logcounts matrix (only keep the 500 most variables log-cpm)
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

colnames(highly_variable_lcpm) = id.infor.3d.ysb$STRAIN

###### exploratory plots #############

pdf("clustering.pdf")
# we apply a hierarchical clustering on the samples and on the rows.
heatmap.2(highly_variable_lcpm, trace="none", main="Top 500 most variable genes across samples")
dev.off()
# splitting by strain

pdf("MDS.pdf")
# Multidimensional scaling (MDS) plot:
plotMDS(y, col=c("red", "blue", "blue", "blue", "red", "red")) # red =3D1
# remember that MDS (and PCA) plots are useful yet partial representations of the data.
dev.off()
# similar variation in both

################ Differential expression analysis
# using highly expressed genes and comparing between heat and control


#Colcounts.CH.group - this file has id, treatment and reads information, you will need something similar. 
# how many reads (library size) did you get for each individual   

# select only highly expressed genes 
head(counts.3d.ysb)
keep <- rowSums(cpm(counts.3d.ysb)>10) >=2 
Hcounts.3d.ysb = counts.3d.ysb[keep,]
dim(Hcounts.3d.ysb) # 8253

# create a DGEList element

ysb3d.str = DGEList(counts=Hcounts.3d.ysb, genes=row.names(Hcounts.3d.ysb), lib.size = id.infor.3d.ysb$total, group = id.infor.3d.ysb$STRAIN)
str(ysb3d.str)
# exploratory plot - Check for outliers. ideally individuals of the same group should group together, but dont worry too much if they done. Good 
plotMDS(ysb3d.str, col=c("red", "blue", "blue", "blue","red", "red"))
# 
# calculate normalisation factors, normalising on library size
ysb3d.str <- calcNormFactors(ysb3d.str)
design <- model.matrix(~id.infor.3d.ysb$STRAIN)

# estimate disperson
ysb3d.str <- estimateDisp(ysb3d.str)
ysb3d.str <- estimateCommonDisp(ysb3d.str)
ysb3d.str <- estimateTagwiseDisp(ysb3d.str)

# testing for differenitally expressed genes
et.ysb3d.str = exactTest(ysb3d.str, pair=levels(ysb3d.str$samples$group))
topTags(et.ysb3d.str) # most differentiated 

dim(et.ysb3d.str) # 8253 3
deg.ysb3d.str <- decideTestsDGE(et.ysb3d.str, adjust.method="BH", p.value=0.05)
summary(deg.ysb3d.str)
#       ST99CH_3D7-ST99CH_3D1
# Down                    2617
# NotSig                  2962
# Up                      2674


sigTags.ysb3d.str = topTags(et.ysb3d.str, n=2617+264)
# write out the list of DEGs

write.table(sigTags.ysb3d.str, file="PATH/out/deg_ysb37dv3d1_high10sig.txt", quote=FALSE, row.names=FALSE)
write.table(et.ysb3d.str, file="PATH/out/eT_ysb3d7v3d1_high10.txt", quote=FALSE, row.names=FALSE)


