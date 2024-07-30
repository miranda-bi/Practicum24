
setwd("C:/Users/Miran/OneDrive/桌面/Practicum/AURC/TWAS_PC1/")

load("TWAS_AURC_v6_PC1_mashr.Rdata")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# library(dplyr)
# library(org.Hs.eg.db)
# library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(snplist)

library(biomaRt)

# Initialize the use of the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Create a list of unique gene identifiers from your dataset
gene_ids <- sub("\\..*$", "", TWAS_v6$gene)
TWAS_v6$gene_ids <- gene_ids

# Retrieve gene location information
gene_locations <- getBM(attributes = c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position'),
                        filters = 'ensembl_gene_id', values = gene_ids, mart = ensembl)

# Merge this location information back into your original dataset
TWAS_v6_res <- merge(TWAS_v6, gene_locations, by.x = 'gene_ids', by.y = 'ensembl_gene_id')

# get the base pair position
TWAS_v6_res$BP <- with(TWAS_v6_res, (start_position + end_position)/2)
TWAS_v6_res$SNP <- paste0("SNP_", seq_len(nrow(TWAS_v6_res))) # create dummy snp
TWAS_v6_res <-na.omit(TWAS_v6_res) # remove nas

library(qqman)

# The Manhattan plot
jpeg("TWAS_Manhattan_Plot.jpeg")
manhattan(TWAS_v6_res, chr="chromosome_name", bp="BP", 
          p="p.value", snp="SNP",
          suggestiveline = -log10(4.17e-06),
          genomewideline = FALSE,
          main="Manhattan Plot for Usual OLS Model for TWAS",
          ylim=c(0, -log10(min(TWAS_v6_res$p.value))+1),
          xlab = "Chromosome",
          ylab = "-log10(P-value)",
          col = c("blue4", "orange3"),
           cex.lab=0.9, cex.main = 0.9)

dev.off()
# change the suggestive sig p-value threshold to 0.05/#genes = 4.17e-06


# QQ plot
jpeg("TWAS_QQ_Plot.jpeg")
qq(TWAS_v6_res$p.value, main = "Q-Q plot for Usual OLS Model for TWAS",
   cex.lab=0.9, cex.main = 0.9)
dev.off()



## LGC

# observed chi-square statistics
obs_chi <- qchisq(TWAS_v6_res$p.value, df = 1, lower.tail = FALSE)

# The median of the observed chi-square statistics
median_obs <- median(obs_chi)

# The median of the expected chi-square distribution with 1 degree of freedom under the null is 0.455
median_exp <- qchisq(0.5, df = 1, lower.tail = FALSE)

# Calculate the genomic control lambda
lgc <- median_obs / median_exp

lgc


### distribution of Estimated Gene Expression Levels
summary(TWAS_v6$Est)
library(ggplot2)

ggplot(TWAS_v6, aes(x = Est)) +
  geom_histogram(binwidth = 5000000, fill = "steelblue", color = "black") +
  labs(title = "Distribution of Estimated Gene Expression Levels",
       x = "Estimated Gene Expression Levels",
       y = "Frequency") +
  theme_minimal()

# remove the outliers


