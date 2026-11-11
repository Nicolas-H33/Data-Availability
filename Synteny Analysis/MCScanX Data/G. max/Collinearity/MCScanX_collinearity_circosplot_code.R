library(circlize)
library(dplyr)
library(stringr)

setwd("/Users/nickh245/Library/CloudStorage/OneDrive-IowaStateUniversity/Objective 1/dNdS analysis/MCScanX_temp")

#Input
gff <- read.table("Glycine_max_filtered.gff", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(gff) <- c("chr", "gene", "start", "end")

collinearity <- readLines("Glycine_max.collinearity")
lines <- collinearity[str_detect(collinearity, "^[0-9]+-")]

#Extract pairs and positions
pairs <- data.frame(
  gene1 = str_extract(lines, "Glyma\\.[0-9]{2}G[0-9]{6}\\.Wm82\\.a6\\.v1"),
  gene2 = str_extract_all(lines, "Glyma\\.[0-9]{2}G[0-9]{6}\\.Wm82\\.a6\\.v1") %>% sapply(`[`, 2),
  evalue = as.numeric(str_extract(lines, "\\S+$")),
  stringsAsFactors = FALSE)

pairs <- pairs %>%
  left_join(gff, by = c("gene1" = "gene")) %>%
  rename(chr1 = chr, start1 = start, end1 = end) %>%
  left_join(gff, by = c("gene2" = "gene")) %>%
  rename(chr2 = chr, start2 = start, end2 = end)

dml_genes <- c("Glyma.10G065900.Wm82.a6.v1", "Glyma.13G151000.Wm82.a6.v1", "Glyma.03G190800.Wm82.a6.v1",
               "Glyma.20G188300.Wm82.a6.v1", "Glyma.10G202150.Wm82.a6.v1")

#Filter
pairs <- pairs %>% filter(!is.na(evalue) & evalue <= 1e-20)

pairs_bg <- pairs %>% filter(!(gene1 %in% dml_genes | gene2 %in% dml_genes))
pairs_dml <- pairs %>% filter(gene1 %in% dml_genes | gene2 %in% dml_genes)

#Chromosomes
chr_sizes <- gff %>%
  group_by(chr) %>%
  summarise(max_end = max(end))

#Circos plot
circos.clear()
set.seed(123)
circos.initialize(factors = chr_sizes$chr, xlim = cbind(0, chr_sizes$max_end))
circos.track(
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    circos.rect(xlim[1], 0, xlim[2], 1, col = "gold2", border = "black", lwd = 1.2)
    circos.text(mean(xlim), 0.5, chr, cex = 0.6, col = "black", facing = "inside", niceFacing = TRUE)
  },
  track.height = 0.08, bg.border = NA)

circos.genomicLink(
  region1 = pairs_bg %>% select(chr1, start1, end1),
  region2 = pairs_bg %>% select(chr2, start2, end2),
  col = adjustcolor("gray65", alpha.f = 0.6),
  border = NA)

circos.genomicLink(
  region1 = pairs_dml %>% select(chr1, start1, end1),
  region2 = pairs_dml %>% select(chr2, start2, end2),
  col = adjustcolor("red", alpha.f = 1),
  border = "red",
  lwd = 3)