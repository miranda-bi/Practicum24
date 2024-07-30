
### This script is intended to be run in Niagara
# salloc -p debug -t 1:00:00 -N 1 --ntasks-per-node=40
# salloc -p compute -t 2:30:00 -N 1 --ntasks-per-node=80

# ml NiaEnv/2019b gcc/9.2.0 openmpi/4.0.3 mkl/2019u4 r/3.6.3
# R

### Define the source directory path and the output directory path
SD="/scratch/o/oespinga/mirandab/GWAS_Logtrans_files/"
OUTDIR <- paste0(Sys.getenv("SCRATCH"),"/GWAS_Logtrans_files/")

PHENO_DIR="/scratch/o/oespinga/reshani0/"

system("ls /scratch/o/oespinga/reshani0")
load(paste0(PHENO_DIR, "demo_v6.Rdata"))
load(paste0(PHENO_DIR, "demo_v8.Rdata"))

# load PCs data
PC_path <- "/project/o/oespinga/oespinga/OAI/Geno/PCs_bdata_n=4113.txt"

# fill = TRUE to handle lines of varying lengths
# strip.white = TRUE to remove extra whitespace
PC_dat <- read.table(PC_path, header = TRUE, sep = "", fill = TRUE,
                     strip.white = TRUE, comment.char = "")

# Check the first few rows of the data
head(PC_dat)

PC1 <- PC_dat[1:3]

# List current loaded data
ls()

### Load phenotype and covariate info 
# prepare data
library(dplyr)

head(demo_v6)
print(unique(demo_v6$V00COHORT))


dat <- demo_v6 %>% 
  filter(V00COHORT!="3: Non-exposed control group") %>% 
  select(ID, V00COHORT, BMI, AGE, SEX, ARUC) %>% 
  rename(AURC = ARUC, FID = ID)

dat$Log_AURC<-log(dat$AURC)
dat<-dat %>% select(-AURC)

head(dat)
print(unique(dat$V00COHORT))

# Merge dat with PC1 by "FID" and keep only matching observations 
merged_dat <- merge(dat, PC1, by = "FID")
head(merged_dat)

# Rearrange the columns
merged_dat <- merged_dat[, c(1, 7:8, 2:6)]

pheno_dat <- merged_dat[, c("FID", "IID", "Log_AURC")] 
covar_dat <- merged_dat[, -8] # Remove outcome "AURC"

# Change categorical columns to numerical
covar_dat <- covar_dat %>%
  mutate(V00COHORT = case_when(
    V00COHORT == '1: Progression' ~ 0,
    V00COHORT == '2: Incidence' ~ 1,)) %>% 
  mutate(SEX = case_when(
    SEX == 'F' ~ 0,
    SEX == 'M' ~ 1,)) 


### Define and save output files

# dt.mat <- model.matrix(~FID+IID+PC1+V00COHORT+BMI+AGE+SEX, data=covar_dat) # change to numeric
# write.table(dt.mat[, -1], file=paste0(OUTDIR, "covar.txt"), col.names=TRUE, row.names = F, quote = F)
# Error when using this method: Error: Fewer tokens than expected on line 2 of --covar file. 
# There is a mismatch in the number of elements between the second line and header line.


# Get covariate data text file
write.table(covar_dat, file=paste0(OUTDIR, "covar.txt"), col.names=TRUE, row.names = F, quote = F)

# Get phenotype data text file
write.table(pheno_dat, file=paste0(OUTDIR, "pheno.txt"), col.names=TRUE, row.names = F, quote = F)



# check if the format is correct 
print(read.table("/scratch/o/oespinga/mirandab/GWAS_Logtrans_files/covar.txt", header = TRUE, nrows = 5))

print(read.table("/scratch/o/oespinga/mirandab/GWAS_Logtrans_files/pheno.txt", header = TRUE, nrows = 5))




##########



