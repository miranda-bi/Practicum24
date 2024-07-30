
# salloc -p debug -t 1:00:00 -N 1 --ntasks-per-node=40
# salloc -p compute -t 2:30:00 -N 1 --ntasks-per-node=80

# ## Copy genotype data files
# cp /project/o/oespinga/oespinga/OAI/Geno/QC_v7/OAI_b37.bim /scratch/o/oespinga/mirandab/GWAS_RobustCov_files/
# cp /project/o/oespinga/oespinga/OAI/Geno/QC_v7/OAI_b37.bed /scratch/o/oespinga/mirandab/GWAS_RobustCov_files/
# cp /project/o/oespinga/oespinga/OAI/Geno/QC_v7/OAI_b37.fam /scratch/o/oespinga/mirandab/GWAS_RobustCov_files/
#   
# ## Copy phenotype and covariate data files
# cp /scratch/o/oespinga/mirandab/Pheno_Cov_files/pheno.txt /scratch/o/oespinga/mirandab/GWAS_RobustCov_files/
# cp /scratch/o/oespinga/mirandab/Pheno_Cov_files/covar.txt /scratch/o/oespinga/mirandab/GWAS_RobustCov_files/

# ml NiaEnv/2019b gcc/9.2.0 openmpi/4.0.3 mkl/2019u4 r/3.6.3
# R

args=(commandArgs(TRUE)) ## retrieve variable from the job call if any
print(args)
if( length(args)>0){
  for(k in 1:length(args)){
    eval(parse(text=args[[k]]))
  }
}
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("snpStats")
# install.packages("nnet")
library(snpStats) # for reading plink files
library(nnet) # for fitting the models

### Define the source directory path and the output directory path
OUTDIR <- paste0(Sys.getenv("SCRATCH"),"/GWAS_RobustCov_files/")

#setwd("C:/Users/Miran/OneDrive/桌面/Practicum/GWAS/RobustCov_files")
setwd("/scratch/o/oespinga/mirandab/GWAS_RobustCov_files/")

gen <- read.plink("OAI_b37.bed", "OAI_b37.bim", "OAI_b37.fam") # split data by CHR
#gen$genotypes: a SnpMatrix with  4113 rows and  2309362 columns
#gen$fam: the sample information, including 'pedigree' (FID) and 'member' (IID), etc
# gen$map: SNP information, including 'chromosome', 'snp.name', 'allele.1&2', etc
phen <- read.table('pheno.txt',header=T)
cov = read.table('covar.txt',header=T)

### merge phenotype and covariate data
dt <- merge(phen, cov)
# dt$V00COHORT: 0 is Progression, 1 is Incidence
# dt$SEX: 0 is female, 1 is male

### Ensure that gen$genotypes and dt datasets have the same subjects and are in the same order
ids.gen <- paste(gen$fam$pedigree,gen$fam$member,sep=".")
ids.dt <- paste(dt$FID,dt$IID,sep=".")

### check that both datasets have the same subjects 
if( !(all(ids.gen%in%ids.dt) & all(ids.dt%in%ids.gen)) ){
  ### to subset the geno info:
  ids_keep <- which(ids.gen %in% ids.dt)
  if( !all.equal(rownames(gen$genotypes),as.character(gen$fam$pedigree)) ) stop("IDs b/t fam and genotype files do not match!")
  
  gen$fam <- gen$fam[ids_keep,] # keep only individuals that are in phenotype data
  gen$genotypes <- gen$genotypes[ids_keep,]# keep only the SNP information for these individuals
} 
# updated 'gen$genotypes' has 1196 rows and 2309362 columns 
# Row name: FID; Col name: snps name

### check if they are in the same order
if( !all(rownames(gen$genotypes)==as.character(dt$FID)) ){
  dt <- dt[order(match(ids.dt, ids.gen)),]
}

# load required packages
#install.packages("lmtest")
library(parallel)
library(doParallel)
library(iterators)
library(foreach)
library(snow)
library(sandwich)
library(lmtest)

cl <- makeCluster(detectCores()-1, type="SOCK", outfile="")
registerDoParallel(cl)

### test the code below by testing for a single column
## i <- sample(NCOL(gen$geno),1); G <- gen$geno[,i]
### or you could also test on a subset of the columns (comment line 58 and comment out lines 56-57)
# indx <- sample(NCOL(gen$geno),500)
# a <- foreach(G=iter(gen$geno[,indx], by='column'), .combine='rbind', .multicombine=T, .inorder=F, .verbose=T, .packages = c("nnet","lmtest","sandwich")) %dopar% {

a <- foreach(G=iter(gen$geno, by='column'), .combine='rbind', .multicombine=T, .maxcombine=1000, .inorder=F, .verbose=F, .packages = c("nnet","lmtest","sandwich")) %dopar% {
  
  snp = colnames(G)
  G=3-as.numeric(G) # reverse the snpStats default coding (according to the number of copies of the MAYOR allele)
  
  ## compute some summary values 
  maf <- mean(G)/2
  N <- length(!is.na(G))
  
  if( maf > 0.05 ){
    ## Fit the alternative model 
    # ordinal linear regression with “usual” OLS standard errors
    # this assumes that the data are uncorrelated and homoscedastic
    md1 <- lm(AURC ~ G + as.numeric(dt$PC1) + as.factor(dt$V00COHORT) + as.numeric(dt$BMI) + 
                      as.numeric(dt$AGE) + as.factor(dt$SEX), data=dt)
    
    # unadjusted standard errors
    se <- summary(md1)$coefficient["G",2]
    
    # t test using the basic “robust” sandwich covariance 
    robust_res <- coeftest(md1, vcov = sandwich)
    
    # Calculate robust standard errors
    robust_se <- robust_res[,2]
    #sqrt(diag(vcovHC(md1, type = "HC0")))
    
    # Fit the null model
    md0 <- lm(AURC ~ 1 + as.numeric(dt$PC1) + as.factor(dt$V00COHORT) + as.numeric(dt$BMI) + 
                as.numeric(dt$AGE) + as.factor(dt$SEX), data=dt)
    
    # LRT: get a p from unadjusted anova
    lrt <- anova(md1, md0)
    
    # Extract coefficients for the SNP effect and use robust SE
    coefs <- coef(md1)["G"]
    se1 <- robust_se["G"]
    t_val <- robust_res["G",3] #coefs/se
    pvals <- robust_res["G",4]
    
    # Assemble results
    res <- unlist(c(maf, coefs, lrt$"Pr(>F)"[2], se, se1, t_val, pvals, N), use.names=F)
    names(res) <- c("maf", "coef", "unadjust_p_value", "unadjust_se", "robust_se", "robust_t", "robust_p_value", "N")
  } else {
    res <- c(maf, rep(NA, 6), N)
    names(res) <- c("maf", "coef", "unadjust_p_value", "unadjust_se", "robust_se", "robust_t", "robust_p_value", "N")
  }
  
  ## return
  data.frame(snp.name=snp, t(res))
} 

stopCluster(cl)

### add  chr, position, allele1 and allele2
a1 = merge(gen$map[,-3],a)
### sort per position
a1 = a1[order(a1$position),]

# remove NAs
a1 <- na.omit(a1)

### save results to text file 
write.table(a1,file=paste0(OUTDIR,"RobustVar_Results"),row.names=F,quote=F,col.names=F)


