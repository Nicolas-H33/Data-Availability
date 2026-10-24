library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)

setwd("/Users/nickh245/Library/CloudStorage/OneDrive-IowaStateUniversity/lab presentations, posters, writing, pictures/My writing")
TF_expression <- read_excel("Supplementary_Tables.xlsx", sheet = "Table S6")

#Parameters
study_prefix <- "fusarium_"
mean         <- 0.25
variance     <- 0.50          
single_study <- FALSE

df <- TF_expression %>%
  filter(grepl(paste0("^", study_prefix), Sample, ignore.case = TRUE))
if (nrow(df) == 0) stop(sprintf("No rows matched prefix '%s'. Check Sample names.", study_prefix))

df <- df %>%
  mutate(
    Susceptibility = factor(Susceptibility),
    `mock/inoculated` = case_when(
      tolower(`mock/inoculated`) == "mock" ~ "mock",
      tolower(`mock/inoculated`) == "inoculated" ~ "inoculated",
      TRUE ~ `mock/inoculated`),
    `mock/inoculated` = factor(`mock/inoculated`, levels = c("mock", "inoculated")),
    Time_h = as.numeric(str_extract(`Time duration`, "\\d+")),
    Time_h = factor(Time_h, levels = sort(unique(Time_h)), ordered = TRUE))

exp_columns <- df %>% select(where(is.numeric)) %>% colnames()

#TFs in columns
base_sel <- c("Susceptibility", "mock/inoculated", "Time_h", exp_columns)
if (single_study) base_sel <- c("Study", base_sel)

long_df <- df %>%
  select(all_of(base_sel)) %>%
  pivot_longer(cols = all_of(exp_columns), names_to = "Gene", values_to = "exp")

#Statistics (mean, variance, & quantiles)
if (single_study) {
  group_means <- long_df %>%
    group_by(Gene, Study, Susceptibility, `mock/inoculated`, Time_h) %>%
    summarize(combo_mean = mean(exp, na.rm = TRUE), .groups = "drop")
} else {
  group_means <- long_df %>%
    group_by(Gene, Susceptibility, `mock/inoculated`, Time_h) %>%
    summarize(combo_mean = mean(exp, na.rm = TRUE), .groups = "drop")
}

gene_stats <- group_means %>%
  group_by(Gene) %>%
  summarize(
    adjusted_mean = mean(combo_mean, na.rm = TRUE),
    between_var   = ifelse(n() > 1, var(combo_mean, na.rm = TRUE), 0),
    n_combos      = n(),
    .groups = "drop")

mean_cut <- quantile(gene_stats$adjusted_mean, probs = mean, na.rm = TRUE)
var_cut  <- quantile(gene_stats$between_var,   probs = variance,  na.rm = TRUE)

filtered_genes <- gene_stats %>%
  mutate(pass_mean = adjusted_mean >= mean_cut, pass_var  = between_var   >= var_cut, pass_both = pass_mean & pass_var)

kept_genes <- filtered_genes %>% filter(pass_both) %>% pull(Gene)
print(kept_genes)

summary <- filtered_genes %>%
  summarize(
    total_genes    = n(),
    pass_mean      = sum(pass_mean),
    pass_var       = sum(pass_var),
    pass_both      = sum(pass_both),
    pass_mean_only = sum(pass_mean & !pass_var),
    pass_var_only  = sum(!pass_mean & pass_var),
    fail_both      = sum(!pass_mean & !pass_var)
  ) %>%
  mutate(pct_pass_both = round(100 * pass_both / total_genes, 1))
print(summary)

#Scatter plot
par(mgp = c(1.8, 0.7, 0))
plot(
  filtered_genes$adjusted_mean, filtered_genes$between_var,
  pch = 19,
  col = ifelse(filtered_genes$pass_both, "black", "grey"),
  xlab = sprintf("Adjusted mean (>= %.0fth pctile = %.3f)", mean*100, mean_cut),
  ylab = sprintf("Between-group variance \n(>= %.0fth pctile = %.3f)", variance*100, var_cut),
  main = "Fusarium: Mean top 75% & Variance top 50%",
  cex.axis = 0.8,
  cex.lab = 0.8,
  cex.main = 0.9)
abline(v = mean_cut, lty = 2, col = "red")
abline(h = var_cut,  lty = 2, col = "red")
legend(
  "top",
  legend = c(paste("Kept =", sum(filtered_genes$pass_both)),
             paste("Dropped =", sum(!filtered_genes$pass_both))),
  pch = 20,
  col = c("black", "grey"),
  bty = "n",
  cex = 0.8)

list(raw_expression = long_df, raw_means = group_means, final = gene_stats,
     variance_quartiles = summary(gene_stats$between_var), mean_quartiles = quantile(gene_stats$adjusted_mean))

#Final table and excel output
kept_TFs <- filtered_genes %>%
  filter(pass_both) %>%                                     
  select(Gene, adjusted_mean, between_var) %>%              
  arrange(desc(between_var))                                

print(kept_TFs, n = Inf)

#write_xlsx(kept_TFs, "Kept_TFs_above_thresholds_fus.xlsx")