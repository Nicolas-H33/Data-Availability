library(gggenomes)
library(dplyr)
library(seqinr)
library(tibble)
library(stringr)

setwd("/Users/nickh245/Library/CloudStorage/OneDrive-IowaStateUniversity/Objective 1/GmROS & GmDME seqs & guides")

#AA seqs
aa <- read.fasta("G.max_5demeth_aa_seqs_updated.fasta")
seqs <- tibble(seq_id = names(aa), length = sapply(aa, length))

#gff3
raw_feats <- read_feats("Introproscan_GmGsAt.gff3")
feats <- raw_feats %>%
  filter(type == "protein_match", grepl("InterPro:IPR", dbxref, ignore.case = TRUE)) %>%
  mutate(
    type = str_extract(dbxref, "InterPro:IPR\\d+"),
    start = as.integer(start),
    end = as.integer(end),
    seq_id = gsub("\\.\\d+$", "", seq_id)) %>%
  select(seq_id, start, end, type) %>%
  filter(!is.na(seq_id), !is.na(start), !is.na(end), !is.na(type),
    seq_id %in% c("GmROS1", "GmROS2", "GmROS3", "GmDME1", "GmDME2"))

feats <- feats %>%
  group_by(seq_id, type) %>%
  filter(!(type == "InterPro:IPR003265" &
             n() > 1 &
             (end - start) > min(end - start))) %>%
  ungroup()

ipr011257_split <- tribble(
  ~seq_id,   ~start, ~end,
  "GmROS1",  1051,   1106,
  "GmROS1",  1461,   1636,
  "GmROS2",  974,    1029,
  "GmROS2",  1433,   1608,
  "GmROS3",  761,    816,
  "GmROS3",  1186,   1361,
  "GmDME1",  825,    881,
  "GmDME1",  1316,   1491,
  "GmDME2",  866,    922,
  "GmDME2",  1357,   1532)

feats <- feats %>%
  filter(!(seq_id %in% ipr011257_split$seq_id &
    type == "InterPro:IPR011257")) %>%
  bind_rows(ipr011257_split %>%
    mutate(type = "InterPro:IPR011257"))

rename <- c(
  "GmROS1" = "GmROS2-1",
  "GmROS2" = "GmROS2-2",
  "GmROS3" = "GmROS1",
  "GmDME1" = "GmDME1-1",
  "GmDME2" = "GmDME1-2")

feats <- feats %>%
  mutate(seq_id = dplyr::recode(seq_id, !!!rename))

#Domain information
domain_labels <- c("InterPro:IPR011257" = "DNA Glycosylase", "InterPro:IPR003265" = "HhH-GPD domain",
                   "InterPro:IPR003651" = "Endonuclease III-like", "InterPro:IPR023170" = "Helix-hairpin-helix",
                   "InterPro:IPR028924" = "zf-CXXC unit", "InterPro:IPR028925" = "Demeter RRM-fold")

domain_colors <- c("DNA Glycosylase" = "blue3", "HhH-GPD domain" = "orange",
                   "Endonuclease III-like" = "cyan", "Helix-hairpin-helix" = "forestgreen",
                   "zf-CXXC unit" = "purple3", "Demeter RRM-fold" = "red2")

feat_list <- lapply(names(domain_labels), function(ipr) {
  feats %>%
  filter(type == ipr) %>%
  mutate(domain = domain_labels[[ipr]])})

names(feat_list) <- domain_labels

#gggenomes plot
gggenomes(seqs = seqs, feats = feat_list) +
  geom_seq(linewidth = 1) +
  geom_seq_label(vjust = 1.2, size = 3.5) +
  geom_feat(data = feats("DNA Glycosylase"),       aes(color = "DNA Glycosylase"), linewidth = 7.5) +
  geom_feat(data = feats("HhH-GPD domain"),        aes(color = "HhH-GPD domain"), linewidth = 5.5) +
  geom_feat(data = feats("Endonuclease III-like"), aes(color = "Endonuclease III-like"), linewidth = 5.5) +
  geom_feat(data = feats("Helix-hairpin-helix"),   aes(color = "Helix-hairpin-helix"), linewidth = 2.8) +
  geom_feat(data = feats("zf-CXXC unit"),          aes(color = "zf-CXXC unit"), linewidth = 5.5) +
  geom_feat(data = feats("Demeter RRM-fold"),      aes(color = "Demeter RRM-fold"), linewidth = 5.5) +
  scale_color_manual(name = "Protein Domains", values = domain_colors) +
  scale_x_continuous(
    breaks = seq(0, 2100, by = 300),
    labels = as.character(seq(0, 2100, by = 300)),
    limits = c(0, 2100),
    name = "Amino Acids",
    expand = c(0, 0)
  ) +
  ggtitle("Protein Domains and Motifs") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 15),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 11),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.text.x  = element_text(size = 10),
    axis.ticks.x = element_line(),
    axis.line.x  = element_line(),
    legend.position = "right",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
    aspect.ratio = 0.6
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 4)))