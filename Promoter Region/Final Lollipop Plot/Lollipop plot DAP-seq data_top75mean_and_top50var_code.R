library(readxl)
library(dplyr)
library(ggplot2)

setwd("~/Library/CloudStorage/OneDrive-IowaStateUniversity/lab presentations, posters, writing, pictures/My writing")

TFs <- read_excel("Supplementary_Tables.xlsx", sheet = "Table S7") %>%
  pull(Gene)

dap_table <- read_excel("Supplementary_Tables.xlsx", sheet = "Table S5")
dap_table <- dap_table %>%
  mutate(
    true_distance_to_tss = case_when(
      gene_strand == "+" ~ peak_center - tss,
      gene_strand == "-" ~ tss - peak_center)) %>%
  filter(true_distance_to_tss >= -1010 & true_distance_to_tss <= 1010)

colnames(dap_table)

order <- c("GmDME1-2", "GmDME1-1", "GmROS1", "GmROS2-2", "GmROS2-1")
dap_table <- dap_table %>%
  filter(gene_name %in% order) %>%
  mutate(gene_name = factor(gene_name, levels = order))

dap_table_filtered <- dap_table %>%
  filter(tf_name %in% TFs)

unique_tf_genes <- dap_table_filtered %>%
  group_by(tf_name) %>%                    
  summarize(num_bind_sites = n(),                   
            num_genes_bound = n_distinct(gene_name),
            bound_genes = paste(unique(gene_name), collapse = ", ")) %>%
  ungroup() %>%
  arrange(desc(num_genes_bound), tf_name)

print(unique_tf_genes)

#TF height proximity
height <- 0.15
increment <- 0.17

dap_table_filtered <- dap_table_filtered %>%
  arrange(gene_name, true_distance_to_tss) %>%
  group_by(gene_name) %>%
  mutate(
    group = cumsum(c(TRUE, diff(true_distance_to_tss) > 20)),
    height = height + (row_number() - match(group, group)) * increment
  ) %>%
  ungroup() %>%
  select(-group)

#ggplot2 plot
ggplot(dap_table_filtered, aes(x = true_distance_to_tss, y = gene_name)) +
  geom_segment(aes(x = -1010, xend = 1010, yend = gene_name), color = "black", linewidth = 2.4) +
  geom_segment(aes(x = true_distance_to_tss, xend = true_distance_to_tss,
                   y = as.numeric(as.factor(gene_name)),
                   yend = as.numeric(as.factor(gene_name)) + height),
               color = "black", linewidth = 1) +
  geom_point(aes(y = as.numeric(as.factor(gene_name)) + height, fill = tf_family),
             shape = 21, color = "black", size = 8.2, stroke = 0.7) +
  scale_x_continuous(
    breaks = c(-1010, -500, 0, 500, 1010),
    labels = c("1kb", "500bp", "TSS", "500bp", "1kb"),
    name = "") +
  scale_y_discrete(expand = expansion(mult = c(0.11, 0.3))) +
  scale_fill_manual(
    name = "TF Families",
    values = c(
      "AP2" = "blue", "B3" = "grey80", "bHLH" = "darkblue",
      "bZIP" = "firebrick", "C2H2" = "lightpink", "Dof" = "magenta1",
      "ERF" = "darkgreen", "G2-like" = "gold2", "LBD" = "lightcyan2",
      "MYB" = "deeppink2", "MYB-related" = "coral1", "SBP" = "oldlace",
      "WRKY" = "orange", "ZF-HD" = "red", "HSF" = "lightseagreen"
    )) +
  ggtitle("TF Binding - Expression Thresholds") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(color = "black", size = 20),
    axis.text.y = element_text(color = "black", size = 15))
