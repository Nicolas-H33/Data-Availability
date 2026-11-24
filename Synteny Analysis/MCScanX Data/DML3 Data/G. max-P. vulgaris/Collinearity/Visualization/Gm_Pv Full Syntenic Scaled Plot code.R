library(dplyr)
library(stringr)
library(ggplot2)

#Input files
setwd("/Users/nickh245/Commonbean_Glycinemax_synteny")

gm <- read.table("Glycine_max_filtered.gff3", sep = "\t")
pv <- read.table("Phaseolus_vulgaris.gff3", sep = "\t")
names(gm) <- names(pv) <- c("chr", "gene", "start", "end")

#Keep main chromosomes
gm <- gm %>% filter(chr %in% sprintf("Gm%02d", 1:20))
pv <- pv %>%
  filter(chr %in% sprintf("Chr%02d", 1:11)) %>%
  mutate(chr = sub("^Chr", "Pv", chr))

#Chromosome sizes
chr_sizes <- function(x, species, order) {
  x %>%
    group_by(chr) %>%
    summarise(length = max(end), .groups = "drop") %>%
    mutate(species = species) %>%
    arrange(match(chr, order))
}

gm_sizes <- chr_sizes(gm, "Gm", sprintf("Gm%02d", 1:20))
pv_sizes <- chr_sizes(pv, "Pv", sprintf("Pv%02d", 1:11))

write.table(bind_rows(pv_sizes, gm_sizes), "Pv_Gm_chr_sizes.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#Collinearity
lines <- readLines("Pv_Gm.collinearity")
headers <- grep("^## Alignment", lines)

parse_block <- function(i) {
  h <- lines[headers[i]]
  if (!grepl("Chr[0-9]{2}&Gm[0-9]{2}|Gm[0-9]{2}&Chr[0-9]{2}", h))
    return(NULL)
  
  last <- if (i < length(headers)) headers[i + 1] - 1 else length(lines)
  x <- lines[(headers[i] + 1):last]
  x <- strsplit(x[grepl("^\\s*[0-9]+-", x)], "\t")
  if (!length(x)) return(NULL)
  
  data.frame(gene1 = sapply(x, `[`, 2), gene2 = sapply(x, `[`, 3),
             orientation = str_extract(h, "(plus|minus)$"))
}

blocks <- Filter(Negate(is.null), lapply(seq_along(headers), parse_block))
pairs <- bind_rows(blocks, .id = "block_id") %>%
  mutate(
    pv_gene = ifelse(grepl("^PHVUL\\.", toupper(gene1)), gene1, gene2),
    gm_gene = ifelse(grepl("^GLYMA\\.", toupper(gene1)), gene1, gene2),
    pv_short = sub("\\.[0-9]+$", "", toupper(pv_gene)),
    gm_short = sub("\\.WM82\\.A6\\.V1$", "", toupper(gm_gene)))

#DMTase gene pairs
pv_dmtase <- toupper(c("Phvul.007g223600", "Phvul.001g186500", "Phvul.003g032300"))
gm_dmtase <- toupper(c("Glyma.10g065900", "Glyma.13g151000", "Glyma.03g190800",
                       "Glyma.20g188300", "Glyma.10g202150"))

dmtase_pairs <- pairs %>%
  filter(pv_short %in% pv_dmtase, gm_short %in% gm_dmtase)

dmtase_plot <- dmtase_pairs %>%
  left_join(pv %>% rename(pv_gene = gene, ref_chr = chr,
                          ref_start = start, ref_end = end), by = "pv_gene") %>%
  left_join(gm %>% rename(gm_gene = gene, target_chr = chr,
                          target_start = start, target_end = end), by = "gm_gene") %>%
  filter(!is.na(ref_chr), !is.na(target_chr)) %>%
  transmute(ref_chr, ref_start = pmax(1, ref_start - 100000),
            ref_end = ref_end + 100000,
            target_chr, target_start = pmax(1, target_start - 100000),
            target_end = target_end + 100000,
            orientation, ref_species = "Pv", target_species = "Gm")

block_region <- function(b) {
  pv_gene <- ifelse(grepl("^PHVUL\\.", toupper(b$gene1)), b$gene1, b$gene2)
  gm_gene <- ifelse(grepl("^GLYMA\\.", toupper(b$gene1)), b$gene1, b$gene2)
  
  a <- pv[match(pv_gene, pv$gene), ]
  g <- gm[match(gm_gene, gm$gene), ]
  
  valid <- !is.na(a$start) & !is.na(g$start)
  a <- a[valid, ]
  g <- g[valid, ]
  
  if (!nrow(a)) return(NULL)
  
  data.frame(ref_chr = a$chr[1], ref_start = min(a$start),
             ref_end = max(a$end), target_chr = g$chr[1],
             target_start = min(g$start), target_end = max(g$end),
             orientation = b$orientation[1],
             ref_species = "Pv", target_species = "Gm")
}

synteny_plot <- bind_rows(lapply(blocks, block_region))

#Scale
gap <- 8000000
gm_width <- sum(gm_sizes$length) + (nrow(gm_sizes) - 1) * gap
scale_pv <- (gm_width - (nrow(pv_sizes) - 1) * gap) / sum(pv_sizes$length)

pv_sizes$length <- round(pv_sizes$length * scale_pv)

scale_regions <- function(x) {
  x$ref_start <- round(x$ref_start * scale_pv)
  x$ref_end <- round(x$ref_end * scale_pv)
  x
}

synteny_plot <- scale_regions(synteny_plot)
dmtase_plot <- scale_regions(dmtase_plot)

#Chromosome tracks
position_chr <- function(x, y) {
  x$xstart <- cumsum(c(1, head(x$length + gap, -1)))
  x$xend <- x$xstart + x$length
  x$y <- y
  x
}

pv_pos <- position_chr(pv_sizes, 0.50)
gm_pos <- position_chr(gm_sizes, 1.75)
positions <- bind_rows(pv_pos, gm_pos)

make_ribbons <- function(data) {
  bind_rows(lapply(seq_len(nrow(data)), function(i) {
    d <- data[i, ]
    a <- positions[positions$chr == d$ref_chr, ]
    g <- positions[positions$chr == d$target_chr, ]
    if (!nrow(a) || !nrow(g)) return(NULL)
    
    rs <- d$ref_start + a$xstart
    re <- d$ref_end + a$xstart
    ts <- d$target_start + g$xstart
    te <- d$target_end + g$xstart
    
    x <- if (d$orientation == "minus")
      c(ts, te, rs, re) else c(ts, rs, re, te)
    y <- if (d$orientation == "minus")
      c(g$y, g$y, a$y + 0.10, a$y + 0.10) else
        c(g$y, a$y + 0.10, a$y + 0.10, g$y)
    
    data.frame(x = x, y = y, group = i)
  }))
}

grey <- make_ribbons(synteny_plot)
red <- make_ribbons(dmtase_plot)

#ggplot2 Plot
p <- ggplot() +
  geom_polygon(data = grey, aes(x, y, group = group), fill = "grey70", colour = NA, alpha = 0.10) +
  geom_polygon(data = red, aes(x, y, group = group), fill = "red", colour = "red", linewidth = 0.4) +
  geom_rect(data = pv_pos, aes(xmin = xstart, xmax = xend, ymin = y, ymax = y + 0.10),
            fill = "purple3", colour = "black", linewidth = 0.4, alpha = 0.85) +
  geom_rect(data = gm_pos, aes(xmin = xstart, xmax = xend, ymin = y, ymax = y + 0.10),
            fill = "gold2", colour = "black", linewidth = 0.4, alpha = 0.85) +
  geom_text(data = pv_pos, aes((xstart + xend)/2, y + 0.20, label = chr),
            size = 2.4, angle = 45, nudge_y = -0.28) +
  geom_text(data = gm_pos, aes((xstart + xend)/2, y + 0.20, label = chr),
            size = 2.4, angle = 45, nudge_y = -0.02) +
  coord_cartesian(ylim = c(0, 2.2), clip = "off") +
  theme_void() +
  theme(legend.position = "none")

print(p)
