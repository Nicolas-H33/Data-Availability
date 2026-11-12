library(dplyr)
library(ggplot2)

setwd("/Users/nickh/Library/CloudStorage/OneDrive-Personal/Nick_MS_Thesis copy/Objective 1/dNdS analysis/MCScanX_temp")
gene_type <- read.table("Glycine_max.gene_type", sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(gene_type) <- c("gene", "dup_type")

#Duplication summary
summary <- gene_type %>%
  group_by(dup_type) %>%
  summarise(count = n()) %>%
  mutate(percentage = round((count / sum(count)) * 100, 2),
    type = factor(dup_type, levels = c(0,1,2,3,4), labels = c("Singleton", "Dispersed", "Proximal", "Tandem", "WGD/Segmental")),
    type_label = paste0(type, " (", count, ")"))

print(summary)

custom_colors <- c("Singleton" = "lightgreen", "Dispersed" = "orange", "Proximal" = "lightpink",
                   "Tandem" = "lightcyan", "WGD/Segmental" = "lavender")

#ggplot2 plot
ggplot(summary, aes(x="", y=count, fill=type)) +
  geom_col(width=1, color="black") +
  coord_polar(theta="y") +
  theme_void() +
  #labs(title="Soybean Genome Duplication Types") +
  scale_fill_manual(
    values = custom_colors,
    labels = summary$type_label) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=10))