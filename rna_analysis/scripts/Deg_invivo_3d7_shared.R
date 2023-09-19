# DEG mapped to 3D7
# following https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html

# need to replace PATH with path to working directory !!

library(edgeR)
library(DESeq2)
library(GO.db) 

library(corrplot)
library(gplots) 

load("deg_invivo_3d7.RData")

counts = read.table("PATH/data/counts_invivo_3d7.txt", header=T)
head(counts)

id.counts = as.data.frame(names(counts))

id.info <- read.table("PATH/data/sample_info_invivo.txt", header=TRUE)
head(id.info)
names(id.info)[1] <- "id"
id.info$id # sample IDs - 

names(id.counts)[1] = "id"
dim(counts) # 11662 48
dim(id.info) # 48 
id.infor.ordered = merge(id.counts, id.info, all.x=TRUE, sort=FALSE)
head(id.infor.ordered)
head(counts)

##### mapping plots - 3D7 percentage mapped is much higher across all the timepoints
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

pdf("raw_counts.pdf")
# plot counts of samples belonging to the same group
plot(counts[,4], counts[,3], main = "raw counts", xlab = "B2", ylab = "B3")
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

#######
head(id.infor.ordered)

id.infor.ordered[,c(1,4)]
names(counts) <- id.infor.ordered$LibraryName

# split 3D versus 1A
head(id.infor.ordered)
id.infor.ordered$order <- seq(1, length(id.infor.ordered$id))
id.infor.3d <- subset(id.infor.ordered, id.infor.ordered$STRAIN=="ST99CH_3D1" | id.infor.ordered$STRAIN=="ST99CH_3D7")
counts.3d <- counts[,c(id.infor.3d$order)]
dim(counts.3d) # 11662 24
getwd()
write.csv(counts.3d, file="counts_invivo_g3d_3d13d7.csv")

id.infor.1a <- subset(id.infor.ordered, id.infor.ordered$STRAIN=="ST99CH_1A5" | id.infor.ordered$STRAIN=="ST99CH_1E4")
counts.1a <- counts[,c(id.infor.1a$order)]
dim(counts.1a) #24
head(counts.1a)
write.csv(counts.1a, file="counts_invivo_g3d_1a5x1e4.csv")

# split into treatments (same time points) and compare between strains
head(counts.3d)
id.infor.3d$order <- seq(1, length(id.infor.3d$id))

id.infor.3d.7d <- subset(id.infor.3d, id.infor.3d$trt=="7d")
counts.3d.7d <- counts.3d[,c(id.infor.3d.7d$order)]
head(counts.3d.7d)

id.infor.3d.12d <- subset(id.infor.3d, id.infor.3d$trt=="12d")
counts.3d.12d <- counts.3d[,c(id.infor.3d.12d$order)]

id.infor.3d.14d <- subset(id.infor.3d, id.infor.3d$trt=="14d")
counts.3d.14d <- counts.3d[,c(id.infor.3d.14d$order)]

id.infor.3d.28d <- subset(id.infor.3d, id.infor.3d$trt=="28d")
counts.3d.28d <- counts.3d[,c(id.infor.3d.28d$order)]


####### 7d ###########
y <- DGEList(counts=counts.3d.7d, group=id.infor.3d.7d$STRAIN)
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

colnames(highly_variable_lcpm) = id.infor.3d.7d$STRAIN

###### exploratory plots #############

# we apply a hierarchical clustering on the samples and on the rows.
heatmap.2(highly_variable_lcpm, trace="none", main="Top 500 most variable genes across samples")

# Multidimensional scaling (MDS) plot:
plotMDS(y, col=c("red", "red", "blue", "blue", "blue", "red"))

# more variation in 3D1, virtually none in 3d1

################ Differential expression analysis
# using highly expressed genes and comparing between heat and control


#Colcounts.CH.group - this file has id, treatment and reads information, you will need something similar. 
# how many reads (library size) did you get for each individual   

# select only highly expressed genes 
head(counts.3d.7d)
head(id.infor.3d.7d)
# rowSums(cpm(counts.3d.7d)>10)
?cpm
keep <- rowSums(cpm(counts.3d.7d)>10) >=2 
Hcounts.3d.7d = counts.3d.7d[keep,]
dim(Hcounts.3d.7d) #8234

# create a DGEList element
?DGEList
g3d.7d.str = DGEList(counts=Hcounts.3d.7d, genes=row.names(Hcounts.3d.7d), lib.size = id.infor.3d.7d$mapped, group = id.infor.3d.7d$STRAIN)
str(g3d.7d.str)
# exploratory plot - Check for outliers. ideally individuals of the same group should group together, but dont worry too much if they done. Good 
plotMDS(g3d.7d.str,  col=c("red", "red", "blue", "blue", "blue", "red"))

# calculate normalisation factors, normalising on library size
#  The default method for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples.
g3d.7d.str <- calcNormFactors(g3d.7d.str) 
design <- model.matrix(~id.infor.3d.7d$STRAIN)
names(g3d.7d.str)
# [1] "counts"  "samples" "genes"  
# estimate disperson

g3d.7d.str <- estimateDisp(g3d.7d.str)
names(g3d.7d.str)
# [1] "counts"             "samples"            "genes"              "common.dispersion"  "trended.dispersion"
# [6] "tagwise.dispersion" "AveLogCPM"          "trend.method"       "prior.df"           "prior.n"           
# [11] "span"  
g3d.7d.str <- estimateCommonDisp(g3d.7d.str)
# names as above + "pseudo.counts"      "pseudo.lib.size"   
g3d.7d.str <- estimateTagwiseDisp(g3d.7d.str) # For routine differential expresion analysis, we use empirical Bayes tagwise dispersions. Note that common dispersion needs to be estimated before estimating tagwise dispersions. 

plotBCV(g3d.7d.str) # coefficient of variation declines with average CPM. This might not be a good model 


# testing for differenitally expressed genes
et.g3d.7d.str = exactTest(g3d.7d.str, pair=levels(g3d.7d.str$samples$group))
topTags(et.g3d.7d.str) # most differentiated 

dim(et.g3d.7d.str) # 8234 3
deg.g3d.7d.str <- decideTestsDGE(et.g3d.7d.str, adjust.method="BH", p.value=0.05)
summary(deg.g3d.7d.str)
#           ST99CH_3D7-ST99CH_3D1
# Down                      67
# NotSig                  7987
# Up                       180

# more differentiated than not!! 
73+141 # 214
sigTags.g3d.7d.str = topTags(et.g3d.7d.str, n=67+180)
# write out the list of DEGs

write.table(sigTags.g3d.7d.str, file="PATH/out/deg_7d_3d7v3d1_high10sig.txt", quote=FALSE, row.names=FALSE)
write.table(et.g3d.7d.str, file="PATH/out/eT_7d_3d7v3d1_high10.txt", quote=FALSE, row.names=FALSE)


####### 12d ###########

y <- DGEList(counts=counts.3d.12d, group=id.infor.3d.12d$STRAIN)
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

colnames(highly_variable_lcpm) = id.infor.3d.12d$STRAIN

###### exploratory plots #############

# we apply a hierarchical clustering on the samples and on the rows.
heatmap.2(highly_variable_lcpm, trace="none", main="Top 500 most variable genes across samples")

# Multidimensional scaling (MDS) plot:
plotMDS(y, col=c("red",  "blue", "blue", "blue","red", "red"))

# more variation in 3D1, virtually none in 3d1


################ Differential expression analysis
# using highly expressed genes and comparing between heat and control


#Colcounts.CH.group - this file has id, treatment and reads information, you will need something similar. 
# how many reads (library size) did you get for each individual   

# select only highly expressed genes 
head(counts.3d.12d)
keep <- rowSums(cpm(counts.3d.12d)>10) >=2 
Hcounts.3d.12d = counts.3d.12d[keep,]
dim(Hcounts.3d.12d) #8473 6

# create a DGEList element

g3d.12d.str = DGEList(counts=Hcounts.3d.12d, genes=row.names(Hcounts.3d.12d), lib.size = id.infor.3d.12d$mapped, group = id.infor.3d.12d$STRAIN)
str(g3d.12d.str)
# exploratory plot - Check for outliers. ideally individuals of the same group should group together, but dont worry too much if they done. Good 
plotMDS(g3d.12d.str,  col=c("red", "blue", "blue", "blue", "red", "red")) # couple of sampels a little different 

# calculate normalisation factors, normalising on library size
g3d.12d.str <- calcNormFactors(g3d.12d.str)
design <- model.matrix(~id.infor.3d.12d$STRAIN)

# estimate disperson
g3d.12d.str <- estimateDisp(g3d.12d.str)
g3d.12d.str <- estimateCommonDisp(g3d.12d.str)
g3d.12d.str <- estimateTagwiseDisp(g3d.12d.str)

# testing for differenitally expressed genes
et.g3d.12d.str = exactTest(g3d.12d.str, pair=levels(g3d.12d.str$samples$group))
topTags(et.g3d.12d.str) # most differentiated 

dim(et.g3d.12d.str) # 8473 3
deg.g3d.12d.str <- decideTestsDGE(et.g3d.12d.str, adjust.method="BH", p.value=0.05)
summary(deg.g3d.12d.str)
#       ST99CH_3D7-ST99CH_3D1
# Down                      21
# NotSig                  8343
# Up                       109

# more differentiated than not!! 

sigTags.g3d.12d.str = topTags(et.g3d.12d.str, n=21+109)
# write out the list of DEGs

write.table(sigTags.g3d.12d.str, file="PATH/out/deg_12d_3d7v3d1_high10sig.txt", quote=FALSE, row.names=FALSE)
write.table(et.g3d.12d.str, file="PATH/out/eT_12d_3d7v3d1_high10.txt", quote=FALSE, row.names=FALSE)



####### 14d ###########


y <- DGEList(counts=counts.3d.14d, group=id.infor.3d.14d$STRAIN)
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

colnames(highly_variable_lcpm) = id.infor.3d.14d$STRAIN

###### exploratory plots #############

# we apply a hierarchical clustering on the samples and on the rows.
heatmap.2(highly_variable_lcpm, trace="none", main="Top 500 most variable genes across samples")

# Multidimensional scaling (MDS) plot:
plotMDS(y, col=c("red", "red", "red", "blue", "blue", "blue"))

################ Differential expression analysis
# using highly expressed genes and comparing between heat and control


#Colcounts.CH.group - this file has id, treatment and reads information, you will need something similar. 
# how many reads (library size) did you get for each individual   

# select only highly expressed genes 
head(counts.3d.14d)
keep <- rowSums(cpm(counts.3d.14d)>10) >=2 
Hcounts.3d.14d = counts.3d.14d[keep,]
dim(Hcounts.3d.14d) #8738 

# create a DGEList element

g3d.14d.str = DGEList(counts=Hcounts.3d.14d, genes=row.names(Hcounts.3d.14d), lib.size = id.infor.3d.14d$mapped, group = id.infor.3d.14d$STRAIN)
str(g3d.14d.str)
# exploratory plot - Check for outliers. ideally individuals of the same group should group together, but dont worry too much if they done. Good 
plotMDS(g3d.14d.str,  col=c("red", "red", "red", "blue", "blue", "blue"))

# calculate normalisation factors, normalising on library size
g3d.14d.str <- calcNormFactors(g3d.14d.str)
design <- model.matrix(~id.infor.3d.14d$STRAIN)

# estimate disperson
g3d.14d.str <- estimateDisp(g3d.14d.str)
g3d.14d.str <- estimateCommonDisp(g3d.14d.str)
g3d.14d.str <- estimateTagwiseDisp(g3d.14d.str)

# testing for differenitally expressed genes
et.g3d.14d.str = exactTest(g3d.14d.str, pair=levels(g3d.14d.str$samples$group))
topTags(et.g3d.14d.str) # most differentiated 

dim(et.g3d.14d.str) # 8738 3
deg.g3d.14d.str <- decideTestsDGE(et.g3d.14d.str, adjust.method="BH", p.value=0.05)
summary(deg.g3d.14d.str)
#       ST99CH_3D7-ST99CH_3D1
# Down                     968
# NotSig                  6782
# Up                       988

# more differentiated than not!! 

sigTags.g3d.14d.str = topTags(et.g3d.14d.str, n=968+988)
# write out the list of DEGs

write.table(sigTags.g3d.14d.str, file="PATH/out/deg_14d_3d7v3d1_high10sig.txt", quote=FALSE, row.names=FALSE)
write.table(et.g3d.14d.str, file="PATH/out/eT_14d_3d7v3d1_high10.txt", quote=FALSE, row.names=FALSE)

####### 28d ###########


y <- DGEList(counts=counts.3d.28d, group=id.infor.3d.28d$STRAIN)
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

colnames(highly_variable_lcpm) = id.infor.3d.28d$STRAIN

###### exploratory plots #############

# we apply a hierarchical clustering on the samples and on the rows.
heatmap.2(highly_variable_lcpm, trace="none", main="Top 500 most variable genes across samples")

# Multidimensional scaling (MDS) plot:
plotMDS(y, col=c("red", "red", "red", "blue", "blue", "blue"))

# more variation in 3D1, virtually none in 3d1

################ Differential expression analysis
# using highly expressed genes and comparing between heat and control


#Colcounts.CH.group - this file has id, treatment and reads information, you will need something similar. 
# how many reads (library size) did you get for each individual   

# select only highly expressed genes 
head(counts.3d.28d)
keep <- rowSums(cpm(counts.3d.28d)>10) >=2 
Hcounts.3d.28d = counts.3d.28d[keep,]
dim(Hcounts.3d.28d) #8738 

# create a DGEList element

g3d.28d.str = DGEList(counts=Hcounts.3d.28d, genes=row.names(Hcounts.3d.28d), lib.size = id.infor.3d.28d$mapped, group = id.infor.3d.28d$STRAIN)
str(g3d.28d.str)
# exploratory plot - Check for outliers. ideally individuals of the same group should group together, but dont worry too much if they done. Good 
plotMDS(g3d.28d.str,  col=c("red", "red", "red", "blue", "blue", "blue"))

# calculate normalisation factors, normalising on library size
g3d.28d.str <- calcNormFactors(g3d.28d.str)
design <- model.matrix(~id.infor.3d.28d$STRAIN)

# estimate disperson
g3d.28d.str <- estimateDisp(g3d.28d.str)
g3d.28d.str <- estimateCommonDisp(g3d.28d.str)
g3d.28d.str <- estimateTagwiseDisp(g3d.28d.str)

# testing for differenitally expressed genes
et.g3d.28d.str = exactTest(g3d.28d.str, pair=levels(g3d.28d.str$samples$group))
topTags(et.g3d.28d.str) # most differentiated 

dim(et.g3d.28d.str) # 7496 3
deg.g3d.28d.str <- decideTestsDGE(et.g3d.28d.str, adjust.method="BH", p.value=0.05)
summary(deg.g3d.28d.str)
#           ST99CH_3D7-ST99CH_3D1
# Down                    1100
# NotSig                  4906
# Up                      1490
# 
# more differentiated than not!! 

sigTags.g3d.28d.str = topTags(et.g3d.28d.str, n=1110+1490)
# write out the list of DEGs

write.table(sigTags.g3d.28d.str, file="PATH/out/deg_28d_3d7v3d1_high10sig.txt", quote=FALSE, row.names=FALSE)
write.table(et.g3d.28d.str, file="PATH/out/eT_28d_3d7v3d1_high10.txt", quote=FALSE, row.names=FALSE)


