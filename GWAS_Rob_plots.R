
setwd("/scratch/o/oespinga/mirandab/GWAS_RobustCov_files/")

read.table(file="RobustVar_Results.Rout")

### Generate summary plots
jpeg("Robust_Manhattan.jpeg")
manhattan(a1, chr="chromosome", bp="position", p="robust_p_value", snp="snp.name",
          col = c("blue4", "orange3"), suggestiveline = -log10(1e-04),
          genomewideline = -log10(4.167e-06), main = "Manhattan Plot with Robust Variance 
          For Linear Regression")

manhattan(a1, chr="chromosome", bp="position", p="unadjust_p_value", snp="snp.name",
          col = c("blue4", "orange3"), suggestiveline = -log10(1e-04),
          genomewideline = -log10(4.167e-06), main = "Manhattan Plot with Robust Variance 
          For Linear Regression")
dev.off()

# Q-Q Plots
jpeg("Robust_QQ.jpeg")
qq(a1$robust_p_value, main = "QQ Plot of P-values For Linear Regression")
qq(a1$unadjust_p_value, main = "QQ Plot of P-values For Linear Regression")
dev.off()

### generate a p-p plot to compare the p values
# Sort the p-values
sorted_robust_p <- sort(a1$robust_p_value, decreasing = FALSE)
sorted_unadjust_p <- sort(a1$unadjust_p_value, decreasing = FALSE)

# Create the P-P plot
jpeg("PP_plot.jpeg")
plot(a1$robust_p_value, a1$unadjust_p_value, xlab = "Robust P-values", ylab = "Unadjusted P-values",
     main = "P-P Plot of Robust vs. Unadjusted P-values")
abline(0, 1, col = "red", lwd = 2)
dev.off()
