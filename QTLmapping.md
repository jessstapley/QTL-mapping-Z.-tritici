# QTL mapping

QTL mapping was perfomred with the R package'qtl2' v2_0.24. I followed the R/qtl2 User Guide (https://kbroman.org/qtl2/assets/vignettes/user_guide.html) and the recomendations in the book "A guide to QTL mapping with R/qtl" (DOI 10.1007/978-0-387-92125-9).

# Data import
First step was to import the data into R. I did this by creating a control file and then using ```read.cross2``` function to read the data into R. I also read in map data as separate object and created a map object (pg.map) that has the maker name, chromosome, genetic position (cM) and physical position (bp).

```
library(qtl2)
write_control_file("data_1a/control.yaml", crosstype = "haploid", geno_file = "genefile.csv",
                   founder_geno_file = NULL, gmap_file = "gmap.csv", pmap_file ="pmap.csv",
                   pheno_file = "phenofile.csv", covar_file = NULL, phenocovar_file = NULL,
                   sex_file = NULL, sex_covar = NULL, sex_codes = NULL,
                   crossinfo_file = NULL, crossinfo_covar = NULL,
                   crossinfo_codes = NULL, geno_codes = NULL, alleles = c("A", "B"),
                   xchr = NULL, sep = ",", na.strings = c("-", "NA"),
                   comment.char = "#", geno_transposed = FALSE,
                   founder_geno_transposed = FALSE, pheno_transposed = FALSE,
                   covar_transposed = FALSE, phenocovar_transposed = FALSE,
                   description = NULL, comments = NULL, overwrite = TRUE)
mapthis <-  read_cross2("data/control.yaml")
pmap <- read.csv("data/pmap.csv")
gmap <- read.csv("data/gmap.csv")
pg.map <-  pmap
pg.map$pos.cM <-gmap$pos
```
Next I inserted pseudomarkers into positions in the grid betweeen real SNPs and then calcuated genotype probabilities.

```
map <- insert_pseudomarkers(mapthis$gmap, step=1)
pr <- calc_genoprob(mapthis, map, error_prob=0.002)
```
# Genome scan and permuation test

Models that take into account the genetic covariance between individuals (i.e. including the kinship matrix in the model) can reduce the false discovery rate in QTL scans and outperform models that do not include this information (https://link.springer.com/article/10.1007/s00122-011-1558-z). So to perform a QTL scan I used a linear mixed model and includded the kinship matrix as a random effect. The first step was to estimate the kinship matrix, the method 'loco' (“leave one chromosome out”) scans each chromosome using a kinship matrix that is calculated using data from all other chromosomes. After running the scan, I perfomed a premutation test to identify the signficance LOD threshold for a QTL peak.

```
kinship_loco <- calc_kinship(pr, "loco", cores=4)
out_pgL <- scan1(pr, pheno=mapthis$pheno[,-1], kinship = kinship_loco, addcovar = NULL, model = "normal")
operm <- scan1perm(pr, pheno=mapthis$pheno[,-1], kinship = kinship_loco, addcovar = NULL, model = "normal", n_perm = 1000, cores=4)

```
# Identify QTL peaks

To identify signficant QTL peaks and thier intervals I looped through each phenotypic trait (ph.list - is a list of all phenotypic traits). The lod threshold (alpha = 0.01) was extracted from permuation results (operm) and then I used the ```find_peaks``` function to identify the QTL exceeding the lod threshold.

```
for (j in ph.list){
    out_pgL <- scan1(pr, pheno=mapthis$pheno[,j], kinship = kinship_loco, addcovar = NULL, model = "normal")
    lodT <- as.data.frame(summary(operm, alpha=0.01))[j]
    lod.pgL <- find_peaks(out_pgL, map, threshold=lodT , drop=1.5)
    }
  
```

# Calculate size of the QTL interval

I calculated the Bayes Credible Interval (95%) around the QTL using 'bayes_int' function to identify the QTL interval. I looped through each peak listed in lod.pgL from previous step. I obtained the peak marker, the physical positions of the marker and the physical size of the interval (bp, Kb) using a map with centimorgan distance and base pair positons for each marker (pg.map).

```
log.pgL.int <- NULL
for (i in 1:length(lod.pgL$lodindex)) {
      ll <- lod.pgL[i,]
      name = ll$lodcolumnll$pheno <- j
      chr.n <-  ll$chr 
      pgmap.chr <-  droplevels(subset(pg.map, pg.map$chr==chr.n))
      pgmap.chr$mk.dist <- abs(ll$pos-pgmap.chr$pos.cM)
      mm <- subset(pgmap.chr, pgmap.chr$mk.dist==min(pgmap.chr$mk.dist))
      ll$marker <- mm$marker[1]
      ll$pos.bp <- min(mm$pos)
      mbb <- subset(pgmap.chr, pgmap.chr$pos.cM>ll$ci_lo & pgmap.chr$pos.cM<ll$ci_hi)
      ll$ci_mkr <- length(mbb$marker)
      ll$ci_lo_bp <-min(mbb$pos)
      ll$ci_hi_bp <- max(mbb$pos)
      ll$ci.interval.kb <- (ll$ci_hi_bp-ll$ci_lo_bp)/1000
      bb95  <-  bayes_int(out_pgL, map, lodcolumn=name, chr=chr.n, prob=0.95)
      ll$b95.ci_lo <- bb95[1,1]
      ll$b95.ci_hi <- bb95[1,3]
      pgmap.chr <-  droplevels(subset(pg.map, pg.map$chr==chr.n))
      mbd95 <- subset(pgmap.chr, pgmap.chr$pos.cM>ll$b95.ci_lo & pgmap.chr$pos.cM<ll$b95.ci_hi)
      ll$bci95_mkr <- length(mbd95$marker)
      ll$bci95_lo_bp <-min(mbd95$pos)
      ll$bci95_hi_bp <- max(mbd95$pos)
      ll$bci95.interval.kb <- (ll$bci95_hi_bp-ll$bci95_lo_bp)/1000
      log.pgL.int <- rbind(log.pgL.int, ll)
      write.csv(log.pgL.int, file=paste0("out/",j,".csv"))
      }
```

# Identify the genes in a QTL interval
Using the physical positions for the QTL intervals I extracted a list of the functional elements in that interval (e.g. genes, exon, cds, etc) from the .gff files. First I read the .gff file into R and then I looped through each QTL peak (and corresponding interval) to identify a list of functional elements. Two tables were created and exported -  one containing just the genes within an interval and the other contained all functional elements in that interval. Also a new file containing each QTL peak, the interval infromation and the number of genes was created (newdat).


```
gff.1a <- read.delim("path/ref.gff3", comment.char = "#", header=FALSE)
names(gff.1a) <- c("chr", "source", "type", "start", "end", "score", "strand", "frame", "attributes")
gff.1a <- droplevels(subset(gff.1a, gff.1a$type!="chromosome")) # remove the line describing chromosome

dat <- read.csv("data/qtl_peak_interval.csv") # file containting all QTL peaks and the 95% bayes confidence intervals (bci)

newdat <- NULL
    for (j in 1:length(dat[,1])){
    dd <- dat[j,]
    gff.chr <- droplevels(subset(gff.1a, gff.1a$chr==dd$chr))
    
    ggb <-  subset(gff.chr, gff.chr$start > dd$bci95_lo_bp)
    ggbb <-  subset(ggb, ggb$end < dd$bci95_hi_bp) 
    genebb <- droplevels(subset(ggbb, ggbb$type=="gene"))
    write.csv(genebb, paste0("out_bci_genes_1a/", dd$trait, "_", dd$marker, "_gene.csv"), row.names=FALSE, quote=FALSE)
    write.csv(ggbb, paste0("out_bci_allann_1a/", dd$trait,  "_", dd$marker, "_all.csv"), row.names=FALSE, quote=FALSE)
    
    dd$n.bci.genes <- length(genebb[,1])
    newdat <- rbind(newdat, dd)
  }

´´´
