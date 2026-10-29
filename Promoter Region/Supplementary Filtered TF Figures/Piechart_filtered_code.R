library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)

#Input files
setwd("~/Library/CloudStorage/OneDrive-IowaStateUniversity/lab presentations, posters, writing, pictures/My writing/Objective 1 Manuscript")
file <- "Supplementary_Tables_Manuscript.xlsx"

#TF threshold
kept_TFs <- read_excel(file, sheet = "Table S7") %>%
  pull(Gene)

#Matched TFs to dMTase genes
tf_summary <- read_excel(file, sheet = "Table S5") %>%
  mutate(
    true_distance_to_tss = if_else(
      gene_strand == "+", peak_center - tss, tss - peak_center),
    tf_family = str_trim(as.character(tf_family))
  ) %>%
  filter(
    between(true_distance_to_tss, -1010, 1010),
    gene_name %in% c("GmDME1-2", "GmDME1-1", "GmROS1",
                     "GmROS2-2", "GmROS2-1"),
    tf_name %in% kept_TFs,
    !is.na(tf_family), tf_family != ""
  ) %>%
  count(tf_family, name = "count", sort = TRUE) %>%
  mutate(
    prop = count / sum(count) * 100,
    tf_family = factor(tf_family, levels = tf_family)
  )

legend_labels <- setNames(
  paste0(tf_summary$tf_family, " (", tf_summary$count, ")"),
  tf_summary$tf_family)

#ggplot2 plot
ggplot(tf_summary, aes(x = "", y = -prop, fill = tf_family)) +
  geom_col(width = 1, color = "black") +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = paste0(round(prop, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 5) +
  scale_fill_manual(
    name = "TF Families",
    values = c(
      "AP2" = "blue", "B3" = "grey80", "bHLH" = "darkblue",
      "bZIP" = "firebrick", "C2H2" = "lightpink", "Dof" = "magenta1",
      "ERF" = "darkgreen", "G2-like" = "gold2", "LBD" = "lightcyan2",
      "MYB" = "deeppink2", "MYB-related" = "coral1", "SBP" = "oldlace",
      "WRKY" = "orange", "ZF-HD" = "red", "HSF" = "lightseagreen"),
    labels = legend_labels) +
  ggtitle("Proportion of TF Families (Filtered for High Expression & Variance)") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 16), legend.position = "right")

print(tf_summary)