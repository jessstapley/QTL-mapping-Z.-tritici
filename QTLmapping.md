# QTL mapping

QTL mapping was perfomred with the R package qtl2 v2_0.24. I followed the R/qtl2 User Guide (https://kbroman.org/qtl2/assets/vignettes/user_guide.html) and the recoomendations in the book "A guide to QTL mapping with R/qtl" (DOI 10.1007/978-0-387-92125-9)

# Data import
First step was to import the data into R. I did this by creating a control file and then using ```read.cross2``` function to read the data into R.

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

```
Next I inserted pseudomarkers into positions in the grid betweeen real SNPs and then calcuated genotype probabilities.

```
map <- insert_pseudomarkers(mapthis$gmap, step=1)
pr <- calc_genoprob(mapthis, map, error_prob=0.002)
```
# Genome scan and permuation test

Models that take into account the genetic covariance between indidivuals (i.e. including the kinship matrix in the model) can reduce the false discovery rate in QTL scans and outperform models that do not include this information in the model (https://link.springer.com/article/10.1007/s00122-011-1558-z). So to perform a QTL scan I used a linear mixed model and includded the kinship matrix as a random effect. The first step was to estimate the kinship matrix, the method 'loco' (“leave one chromosome out”) scans each chromosome using a kinship matrix that is calculated using data from all other chromosomes. After running the scan, I perfomed a premutation test to identify the signficance LOD threshold for a QTL peak.

```
kinship_loco <- calc_kinship(pr, "loco", cores=4)
out_pgL <- scan1(pr, pheno=mapthis$pheno[,-1], kinship = kinship_loco, addcovar = NULL, model = "normal")
operm <- scan1perm(pr, pheno=mapthis$pheno[,-1], kinship = kinship_loco, addcovar = NULL, model = "normal", n_perm = 1000, cores=4)

```
# Identify QTL peaks and intervals

To identify signficance QTL peaks and thier intervals I looped through each phenotypic trait (ph.list - is a list of all phenotypic traits). The lod threshold (alpha = 0.01) was extracted from permuation results (operm) and then I used the ```find_peaks``` function to identify the QTL exceeding the lod threshold.

```
for (j in ph.list){
    out_pgL <- scan1(pr, pheno=mapthis$pheno[,j], kinship = kinship_loco, addcovar = NULL, model = "normal")
    lodT <- as.data.frame(summary(operm, alpha=0.01))[j]
    
  lod.pgL <- find_peaks(out_pgL, map, threshold=lodT , drop=1.5)
  }
  
```
