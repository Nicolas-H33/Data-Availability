library(gggenomes)
library(Biostrings)
library(universalmotif)
library(dplyr)
library(tibble)
library(grid)

setwd("~/Library/CloudStorage/OneDrive-IowaStateUniversity/Objective 1/GmROS & GmDME seqs & guides/Soybean fasta and alignment files/Gm_Meme motif Data")
fasta    <- "G.max_5demeth_aa_seqs copy.fasta"
meme_txt <- "meme.txt"

aa <- readAAStringSet(fasta)
names(aa) <- sub(" .*", "", names(aa))

rename <- c("GmROS2-1", "GmROS2-2", "GmROS1", "GmDME1-1", "GmDME1-2")
stopifnot(length(names(aa)) == length(rename))
names(aa) <- rename
Newnames <- rename

seqs <- tibble(seq_id = names(aa), length = width(aa)) |>
  mutate(seq_id = factor(seq_id, levels = Newnames))

#All hits above 1e-10
motifs <- read_meme(meme_txt)
hits <- scan_sequences(
  motifs, aa,
  threshold = 1e-10, threshold.type = "pvalue",
  RC = FALSE)

hits_table <- as.data.frame(hits)

feat_tbl <- hits_table |>
  transmute(
    seq_id = sequence,
    start,
    end   = stop,
    name  = motif,
    pval  = pvalue
  ) |>
  arrange(seq_id, start)

#Motif labels
motif_levels <- vapply(motifs, function(m) m@name, character(1))

feat_tbl <- feat_tbl |>
  mutate(
  motif_idx = match(name, motif_levels),
  label = factor(paste0("Motif ", motif_idx),
                        levels = paste0("Motif ", 1:20)))

motif_counts <- feat_tbl |>
  group_by(seq_id) |>
  summarise(motif_count = n(), .groups = "drop")
print(motif_counts)

colors <- c(
  "Motif 1"  = "red", "Motif 2"  = "blue", "Motif 3"  = "yellow",
  "Motif 4"  = "lightblue", "Motif 5"  = "forestgreen", "Motif 6"  = "cyan",
  "Motif 7"  = "pink", "Motif 8"  = "purple", "Motif 9"  = "grey",
  "Motif 10" = "gold", "Motif 11" = "orange", "Motif 12" = "lavender",
  "Motif 13" = "tan", "Motif 14" = "magenta", "Motif 15" = "palevioletred",
  "Motif 16" = "palegreen", "Motif 17" = "wheat", "Motif 18" = "mediumpurple",
  "Motif 19" = "darkcyan", "Motif 20" = "darkred")

feat_tbl$seq_id <- factor(feat_tbl$seq_id, levels = Newnames)

#gggenomes plot
p <- gggenomes(seqs = seqs, feats = feat_tbl) +
  geom_seq(linewidth = 1) +
  geom_seq_label(vjust = 1.2, size = 3.5) +
  geom_feat(aes(color = label), size = 6) +
  scale_color_manual(
    name   = "Motif",
    values = colors,
    limits = paste0("Motif ", 1:20),
    drop   = FALSE
  ) +
  scale_x_continuous(
    breaks = seq(0, 2100, by = 300),
    labels = as.character(seq(0, 2100, by = 300)),
    limits = c(0, 2100),
    name   = "Amino Acids",
    expand = c(0, 0)
  ) +
  guides(color = guide_legend(
    override.aes = list(linewidth = 4.3),
    keyheight = unit(4, "mm"),
    keywidth  = unit(8, "mm"),
    ncol = 1
  )) +
  ggtitle("Protein Domains and Motifs") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 15),
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 10),
    legend.spacing.y = unit(3, "mm"),
    legend.spacing.x = unit(2, "mm"),
    legend.key.height= unit(4, "mm"),
    legend.key.width = unit(12, "mm"),
    legend.position  = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x      = element_line(),
    axis.ticks.x     = element_line(color = "black"),
    axis.title.x     = element_text(size = 13),
    axis.text.x      = element_text(size = 11, color = "black"),
    axis.title.y     = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    plot.margin      = margin(t = 10, r = 40, b = 10, l = 10),
    aspect.ratio     = 0.6)

print(p)
