library(stringr)
library(ggplot2)
library(gggenomes)
library(tibble)

#Input files
base <- "/Users/nickh245/Genomes and Genomic Data_GmGsAt"
at_file <- paste0(base, "/Full Arabidopsis genome TAIR10/Athaliana_167_TAIR10.gene.gff3")
pv_file <- paste0(base, "/common bean partial/Pvulgaris_442_v2.1.gene.gff3")
gm_file <- paste0(base, "/Full Soybean genome Wm82.a6.v1/Soybean.Wm82.a6.v1.gene.gff3/Phytozome/PhytozomeV13/Gmax/Wm82.a6.v1/annotation/Gmax_880_Wm82.a6.v1.gene.gff3")

#Collinear gene pairs
anchors <- tribble(
  ~comparison, ~ref, ~soy,
  "At_Gm07","AT4G33990","Glyma.07G211900",
  "At_Gm07","AT4G34000","Glyma.07G213100",
  "At_Gm07","AT4G34020","Glyma.07G213200",
  "At_Gm07","AT4G34030","Glyma.07G213500",
  "At_Gm07","AT4G34040","Glyma.07G214600",
  "At_Gm07","AT4G34050","Glyma.07G214700",
  "At_Gm07","AT4G34100","Glyma.07G215200",
  "At_Gm07","AT4G34120","Glyma.07G215400",
  "At_Gm07","AT4G34131","Glyma.07G215500",
  "At_Gm02","AT4G33990","Glyma.02G102000",
  "At_Gm02","AT4G34040","Glyma.02G103300",
  "At_Gm02","AT4G34050","Glyma.02G103500",
  "At_Gm02","AT4G34100","Glyma.02G103800",
  "At_Gm02","AT4G34110","Glyma.02G103900",
  "At_Gm02","AT4G34120","Glyma.02G104200",
  "At_Gm02","AT4G34131","Glyma.02G104300",
  "Pv_Gm09","Phvul.003G032100","Glyma.09G118300",
  "Pv_Gm09","Phvul.003G032500","Glyma.09G117900",
  "Pv_Gm09","Phvul.003G033100","Glyma.09G117600",
  "Pv_Gm09","Phvul.003G033200","Glyma.09G117400",
  "Pv_Gm09","Phvul.003G033300","Glyma.09G117300",
  "Pv_Gm09","Phvul.003G033400","Glyma.09G117200")

at02 <- subset(anchors, comparison == "At_Gm02")
at07 <- subset(anchors, comparison == "At_Gm07")
pv09 <- subset(anchors, comparison == "Pv_Gm09")

read_genes <- function(file) {
  x <- read.delim(file, header = FALSE, comment.char = "#", quote = "")
  names(x) <- c("chr","source","type","start","end","score","strand","phase","attr")
  x <- x[x$type == "gene", ]
  
  name <- str_match(x$attr, "Name=([^;]+)")[,2]
  id <- str_match(x$attr, "ID=([^;]+)")[,2]
  x$gene <- ifelse(!is.na(name), name, id)
  x$gene <- sub("^gene:", "", x$gene)
  x$gene <- sub("\\.(Wm82\\.a6\\.v1|v2\\.1|TAIR10)$", "", x$gene)
  
  x[,c("chr","gene","start","end","strand")]
}

at <- read_genes(at_file)
pv <- read_genes(pv_file)
gm <- read_genes(gm_file)

#Genomic regions
get_region <- function(df, genes, pad) {
  hits <- df[df$gene %in% genes, ]
  if (!nrow(hits)) stop("Requested genes not found")
  
  chr <- names(sort(table(hits$chr), decreasing = TRUE))[1]
  hits <- hits[hits$chr == chr, ]
  lo <- max(1, min(hits$start) - pad)
  hi <- max(hits$end) + pad
  
  df[df$chr == chr & df$end >= lo & df$start <= hi, ]
}

regions <- list(
  "Soybean Gm02" = get_region(gm, at02$soy, 15000),
  "Arabidopsis Chr4" = get_region(at, unique(c("AT4G34060",at02$ref,at07$ref)), 10000),
  "Soybean Gm07" = get_region(gm, at07$soy, 15000),
  "Common bean Chr03" = get_region(pv, unique(c("Phvul.003G032300",pv09$ref)), 10000),
  "Soybean Gm09" = get_region(gm, pv09$soy, 15000))

focal <- c("AT4G34060", "Phvul.003G032300")
syntenic <- unique(c(anchors$ref, anchors$soy))

#Chromosomes
prep <- function(df, name) {
  start <- min(df$start)
  span <- max(df$end) - start
  
  data.frame(seq_id = name, feat_id = df$gene,
             start = (df$start - start) / span * 100000 + 1,
             end = (df$end - start) / span * 100000 + 1,
             strand = df$strand, type = "CDS",
             gene_class = ifelse(df$gene %in% focal, "DML3",
                                 ifelse(df$gene %in% syntenic, "Syntenic", "Other")))
}

genes_all <- do.call(rbind, Map(prep, regions, names(regions)))

#Relative gene positions
rel_pos_map <- c(
  "AT4G33990"="-6", "AT4G34000"="-5", "AT4G34020"="-4",
  "AT4G34030"="-3", "AT4G34040"="-2", "AT4G34050"="-1",
  "AT4G34100"="+4", "AT4G34110"="+5", "AT4G34120"="+6",
  "AT4G34131"="+7",
  "Glyma.07G211900"="-6", "Glyma.07G213100"="-5",
  "Glyma.07G213200"="-4", "Glyma.07G213500"="-3",
  "Glyma.07G214600"="-2", "Glyma.07G214700"="-1",
  "Glyma.07G215200"="+4", "Glyma.07G215400"="+6",
  "Glyma.07G215500"="+7",
  "Glyma.02G102000"="-6", "Glyma.02G103300"="-2",
  "Glyma.02G103500"="-1", "Glyma.02G103800"="+4",
  "Glyma.02G103900"="+5", "Glyma.02G104200"="+6",
  "Glyma.02G104300"="+7",
  "Phvul.003G032100"="-2", "Phvul.003G032500"="+2",
  "Phvul.003G033100"="+7", "Phvul.003G033200"="+8",
  "Phvul.003G033300"="+9", "Phvul.003G033400"="+10",
  "Glyma.09G118300"="-2", "Glyma.09G117900"="+2",
  "Glyma.09G117600"="+7", "Glyma.09G117400"="+8",
  "Glyma.09G117300"="+9", "Glyma.09G117200"="+10")

genes_all$relative_position <- unname(rel_pos_map[genes_all$feat_id])

seqs_all <- data.frame(
  seq_id = names(regions),
  length = sapply(names(regions), function(n)
    max(genes_all$end[genes_all$seq_id == n])))

#Construct exact syntenic links
make_links <- function(pairs, ref, soy, ref_name, soy_name, label) {
  r <- ref[match(pairs$ref, ref$feat_id), ]
  s <- soy[match(pairs$soy, soy$feat_id), ]
  
  if (anyNA(r$feat_id) || anyNA(s$feat_id))
    stop(paste("Missing validated anchor:", label))
  
  data.frame(seq_id = ref_name, start = r$start, end = r$end, seq_id2 = soy_name, start2 = s$start, end2 = s$end, comparison = label)
}

links_all <- rbind(
  make_links(at02, subset(genes_all, seq_id == "Arabidopsis Chr4"),
             subset(genes_all, seq_id == "Soybean Gm02"), "Arabidopsis Chr4", "Soybean Gm02", "At-Gm02"),
  make_links(at07, subset(genes_all, seq_id == "Arabidopsis Chr4"),
             subset(genes_all, seq_id == "Soybean Gm07"), "Arabidopsis Chr4", "Soybean Gm07", "At-Gm07"),
  make_links(pv09, subset(genes_all, seq_id == "Common bean Chr03"),
             subset(genes_all, seq_id == "Soybean Gm09"), "Common bean Chr03", "Soybean Gm09", "Pv-Gm09"))

#Track labels
name_map <- c("Soybean Gm02"="Gm02 (#4385)", "Arabidopsis Chr4"="At04 (#4385, #4590)", "Soybean Gm07"="Gm07 (#4590)",
              "Common bean Chr03"="Pv03 (#2572)", "Soybean Gm09"="Gm09 (#2572)")

seqs_all$seq_id <- unname(name_map[seqs_all$seq_id])
genes_all$seq_id <- unname(name_map[genes_all$seq_id])
links_all$seq_id <- unname(name_map[links_all$seq_id])
links_all$seq_id2 <- unname(name_map[links_all$seq_id2])

#gggenomes plot
p <- gggenomes(seqs = seqs_all, genes = genes_all, links = links_all) +
  geom_link(aes(fill = comparison), alpha = 0.30, color = NA) +
  geom_seq(color = "grey30", size = 0.8) +
  geom_seq_label(aes(label = seq_id), size = 3.5) +
  geom_gene(aes(fill = gene_class), color = "grey25", size = 5) +
  geom_gene_tag(aes(label = ifelse(feat_id == "AT4G34060", "AtDML3",
                                   ifelse(feat_id == "Phvul.003G032300",
                                          "PvDML3", ""))), size = 3, nudge_y = 0.15) +
  geom_gene_tag(aes(label = ifelse(gene_class == "Syntenic", relative_position, "")), size = 3, nudge_y = 0.15) +
  scale_fill_manual(
    values = c("DML3"="#D73027", "Syntenic"="#4575B4", "Other"="grey80", 
               "At-Gm02"="#FDAE61", "At-Gm07"="#74ADD1", "Pv-Gm09"="#66BD90"),
    breaks = c("DML3","Syntenic","Other","At-Gm02","At-Gm07","Pv-Gm09"),
    labels = c("DML3","Syntenic Genes","Other","At-Gm02","At-Gm07","Pv-Gm09")) +
  labs(title = "DML3 Microsynteny Scaled", fill = NULL) +
  theme_void() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", vjust = -4,
                                  hjust = 0.5, size = 20))
print(p)