### This script is intended to be run in Niagara
# salloc -p debug -t 1:00:00 -N 1 --ntasks-per-node=40
# salloc -p compute -t 2:30:00 -N 1 --ntasks-per-node=80
# ml NiaEnv/2019b gcc/9.2.0 openmpi/4.0.3 mkl/2019u4 r/3.6.3
# R

### Define the source directory path and the output directory path
SD="/gpfs/fs0/project/o/oespinga/oespinga/OAI/Geno/TWAS/"
SD2="/scratch/o/oespinga/mirandab/Pheno_Cov_files/"
OUTDIR <- paste0(Sys.getenv("SCRATCH"),"/TWAS_files/")


phen <- read.table(paste0(SD2, 'pheno.txt'),header=T)
phen$AURC <- log(phen$AURC)
cov = read.table(paste0(SD2, 'covar.txt'),header=T)

dat <- merge(phen, cov, by="FID")

# Extract outcome variable name and covariate names
out <- names(dat)[3]
covs <- colnames(dat)[5:9]

library(data.table)
library(aod)
library(dplyr)

# Define folder
folder <- "mashr"

# Set the working directory 
setwd(paste0(SD,folder))
ptm <- Sys.time() 

### run TWAS in parallel for each row
Log_TWAS_v6 <- parallel::mclapply(1:22, function(chr){
  # chr=1
  predfile <- grep(paste0("predict_Whole_Blood_chr",chr,".txt"),dir(), value = T)
  predxcanchr <- read.table(predfile, header = T) # actual predicted gene expression
  
  predxcanchr$FID <- as.numeric(sub("_.*", "", predxcanchr$FID))
  
  # Merge the phenotype/covariate data with the prediction data 
  dat.twas.chr <- merge(dat, predxcanchr, by="FID")
  dat.twas.chr <- na.omit(dat.twas.chr) # Remove rows with missing values
  
  # Get gene names from the prediction data
  genenames <- names(predxcanchr)[-(1:2)]
  
  # Perform regression analysis for each gene
  #system.time({
  chrres <- suppressWarnings(as.data.frame(rbindlist( lapply(genenames, function(geneid){
    ## geneid <- sample(genenames,1)
    # cat(geneid,"\n")
    # Construct the formula for the regression model with the outcome, the gene expression, and covariates
    form <- as.formula(paste0(out," ~ ", geneid," + ",paste(covs, collapse = "+")))
    # Create a model frame based on the formula and the data
    mod <- model.frame(form, data=dat.twas.chr)
    y = mod[,out] # Extract the outcome variable from the model frame
    x = model.matrix(form, mod) # Create a model matrix for the predictors
    
    # fit linear regression if there is variation in the gene expression predictor
    if( var(mod[,geneid])>0 ){
      fit <- tryCatch(lm(form, data=dat.twas.chr), error=function(e){print(e); NULL})
    }else fit <- NULL
    
    # extract results
    if( is.null(fit) ){
      # If the model did not fit, create a placeholder data frame with NA values
      fitres <- data.frame(gene=geneid, term=NA, NA, NA, NA, NA)
    }else{
      # Get the summary of the fitted model coefficients
      fitres <- summary(fit)$coef
      # Identify rows in the summary that correspond to the gene expression predictor
      idrows <- grepl(geneid,rownames(fitres))
      # If there are two coefficients for the gene (expression and interaction), perform a Wald test
      if( sum(idrows)==2 ){
        p2df <- tryCatch(wald.test(b=coef(fit), Sigma=vcov(fit), Terms=which(idrows))$result$chi2, error = function(e){print(e); NA})
        # includes columns for "gene" and "term",
        # and specific rows from the original fitres data frame where the condition specified by idrows is satisfied
        fitres <- data.frame(gene=geneid, term=c("expr","int"), fitres[idrows,,drop=FALSE])
        fitres <- rbind(fitres, data.frame(gene=geneid, term="2df",Estimate=NA, Std..Error=NA, z.value=p2df[1], Pr...z..=p2df[3]))
        # If there is only one coefficient for the gene, label it accordingly
      }else if( sum(idrows)==1 ) {
        fitres <- data.frame(gene=geneid, term=ifelse(rownames(fitres)[idrows]==geneid,"expr","int"), fitres[idrows,,drop=FALSE])
        # If something unexpected happens, stop the execution and print an error message
      }else stop("Something unusual happened in model fit, gene = ",geneid)
      # Remove row names from the results data frame
      rownames(fitres) <- NULL
    }
    # Rename the columns of the results data frame
    names(fitres)[-c(1:2)] <- c("Est","StdErr","Statistic","p.value")
    
    # return the results data frame for the current gene
    return(fitres)
  }))))
  #})
  
  # set results as a data frame
  data.frame(chrres)
  
}, mc.cores=parallel::detectCores() )

# Combine all TWAS results into a single data frame
Log_TWAS_v6 <- as.data.frame(rbindlist(Log_TWAS_v6))

save(Log_TWAS_v6, file=paste0(OUTDIR,"Log_TWAS_AURC_v6_PC1_",folder,".RData"))

