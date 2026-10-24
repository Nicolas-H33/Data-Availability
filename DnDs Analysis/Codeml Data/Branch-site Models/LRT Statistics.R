#Log-likelihood values from the nulls and alternative hypothesis in the Codeml runs
branch_site <- data.frame(
  gene = c("GmROS2-1", "GmROS2-2", "GmROS1", "GmDME1-1", "GmDME1-2"),
  L_null = c(-15742.43, -15747.17, -15737.6, -15744.11, -15743.84),
  L_alt  = c(-15742.43, -15747.17, -15735.25, -15744.11, -15743.32))

#LRT statistic and p-value
branch_site$LRT <- 2 * (branch_site$L_alt - branch_site$L_null)
branch_site$p_value <- 1 - pchisq(branch_site$LRT, df = 1)

#Significance (0.05 threshold)
branch_site$significant <- ifelse(branch_site$p_value < 0.05, "Yes", "No")

print(branch_site)
