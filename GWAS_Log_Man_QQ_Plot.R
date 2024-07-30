
### This script is intended to be run in Niagara
# salloc -p debug -t 1:00:00 -N 1 --ntasks-per-node=40
# salloc -p compute -t 2:30:00 -N 1 --ntasks-per-node=80

# ml NiaEnv/2019b gcc/9.2.0 openmpi/4.0.3 mkl/2019u4 r/3.6.3
# R

SD <- "/scratch/o/oespinga/mirandab/GWAS_Logtrans_files/"
setwd("/gpfs/fs0/scratch/o/oespinga/mirandab/Results/")

results_lin <- read.table(paste0(SD, "linear_results.assoc_2.linear"), head=TRUE)

# Manhattan plot and Q-Q plot
#install.packages("qqman")
library("qqman")
jpeg("Log_GWAS_Manhattan.jpeg")
manhattan(results_lin,chr="CHR",bp="BP",p="P",snp="SNP",
          suggestiveline = -log10(1e-05),
          genomewideline = FALSE, 
          main = "Manhattan plot for Log-linear Model for GWAS",
          col = c("steelblue4", "lightskyblue3"))
dev.off()

jpeg("Log_GWAS_QQ-Plot.jpeg")
qq(results_lin$P, main = "Q-Q plot for Log-linear Model for GWAS")
dev.off()



## LGC

# observed chi-square statistics
obs_chi <- qchisq(results_lin$P, df = 1, lower.tail = FALSE)

# The median of the observed chi-square statistics
median_obs <- median(obs_chi)

# The median of the expected chi-square distribution with 1 degree of freedom under the null is 0.455
median_exp <- qchisq(0.5, df = 1, lower.tail = FALSE)

# Calculate the genomic control lambda
lgc <- median_obs / median_exp






# Association
results_ass <- read.table(paste0(SD, "assoc_2.results.qassoc"), head=TRUE)

jpeg("Log_Assoc_Manhattan.jpeg")
manhattan(results_ass,chr="CHR",bp="BP",p="P",snp="SNP",
          main = "Manhattan plot: log association")
dev.off()  

jpeg("Log_Assoc_QQ-Plot.jpeg")
qq(results_ass$P, main = "Q-Q plot: log association")
dev.off()

