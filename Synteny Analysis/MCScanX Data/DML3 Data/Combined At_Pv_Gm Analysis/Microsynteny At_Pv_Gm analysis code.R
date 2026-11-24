#Input files
at <- "/Users/nickh245/Arabidopsis_Glycinemax_synteny"
pv <- "/Users/nickh245/Commonbean_Glycinemax_synteny"

files <- list(
  at1 = c(paste0(at, "/At_Gm.collinearity"), paste0(at, "/At_Gm.gff")),
  at3 = c(paste0(at, "/Arabidopsis_Glycinemax_synteny2/At_Gm2.collinearity"),
          paste0(at, "/Arabidopsis_Glycinemax_synteny2/At_Gm2.gff")),
  pv = c(paste0(pv, "/Pv_Gm.collinearity"), paste0(pv, "/Pv_Gm.gff")))

targets <- c(at1 = "AT4G34060", at3 = "AT4G34060",
             pv = "Phvul.003G032300")

data <- lapply(files, function(f) list(
  col = readLines(f[1]),
  gff = read.table(f[2], col.names = c("chr", "gene", "start", "end"))))

#Extract 10 neighboring genes on each side of DML3
get_neighbors <- function(gff, target) {
  i <- match(target, gff$gene)
  if (is.na(i)) stop(paste(target, "not found"))
  
  x <- gff[gff$chr == gff$chr[i], ]
  x <- x[order(x$start, x$end), ]
  i <- match(target, x$gene)
  x <- x[max(1, i-10):min(nrow(x), i+10), ]
  x$relative_position <- seq_len(nrow(x)) - match(target, x$gene)
  x
}

#Identify soybean matches
find_matches <- function(neighbors, lines) {
  out <- list()
  
  for (i in seq_len(nrow(neighbors))) {
    hits <- grep(neighbors$gene[i], lines, value = TRUE)
    hits <- hits[grepl("Glyma\\.", hits)]
    
    for (hit in hits) {
      soy <- strsplit(trimws(hit), "\\s+")[[1]]
      soy <- soy[grepl("^Glyma\\.", soy)]
      if (!length(soy)) next
      out[[length(out) + 1]] <- data.frame(
        reference_gene = neighbors$gene[i],
        relative_position = neighbors$relative_position[i],
        soybean_gene = soy[1],
        soybean_chr = sub("^Glyma\\.([0-9]+)G.*", "Gm\\1", soy[1]),
        block = sub("^\\s*([0-9]+)-.*", "\\1", hit))
    }}
  if (!length(out)) return(data.frame())
  do.call(rbind, out)
}

matches <- lapply(names(data), function(n)
  find_matches(get_neighbors(data[[n]]$gff, targets[n]), data[[n]]$col))
names(matches) <- names(data)

#Summary of candidate blocks 
match_summary <- lapply(matches, function(x) {
  if (nrow(x) == 0) return(data.frame())
  aggregate(reference_gene ~ block + soybean_chr, data = x,
            FUN = function(g) length(unique(g)))
})

print(match_summary)

#Manually input candidates
blocks <- data.frame(
  dataset = c("at1", "at1", "at3", "at3", "pv"),
  species = c(rep("Arabidopsis thaliana", 4), "Phaseolus vulgaris"),
  run = c(".1", ".1", ".3", ".3", NA),
  alignment = c("4590", "4385", "4590", "4385", "2572"),
  chromosome = c("Gm07", "Gm02", "Gm07", "Gm02", "Gm09"),
  comparison = c("Arabidopsis_.1_vs_Gm07", "Arabidopsis_.1_vs_Gm02",
                 "Arabidopsis_.3_vs_Gm07", "Arabidopsis_.3_vs_Gm02",
                 "CommonBean_vs_Gm09"))

#Summary and anchor pairs tables
results <- lapply(seq_len(nrow(blocks)), function(i) {
  b <- blocks[i, ]
  x <- matches[[b$dataset]]
  x <- x[x$block == b$alignment, ]
  
  up <- unique(x$reference_gene[x$relative_position < 0])
  down <- unique(x$reference_gene[x$relative_position > 0])
  
  summary <- data.frame(Reference_species = b$species, Arabidopsis_protein_run = b$run,
    Reference_DML3 = targets[b$dataset], MCScanX_alignment = b$alignment,
    Soybean_chromosome = b$chromosome, Flanking_genes_matched = length(unique(x$reference_gene)),
    Upstream_genes_matched = length(up), Downstream_genes_matched = length(down),
    Genes_on_both_sides = ifelse(length(up) > 0 & length(down) > 0, "YES", "NO"),
    DML3_MCScanX_anchor = ifelse(any(grepl(targets[b$dataset], data[[b$dataset]]$col)), "YES", "NO"))
  
  anchors <- data.frame(comparison = b$comparison, protein_run = b$run,
    alignment = b$alignment, x[, c("reference_gene", "relative_position", "soybean_gene", "soybean_chr")])
  
  list(summary = summary, anchors = anchors)
})

summary_table <- do.call(rbind, lapply(results, `[[`, "summary"))
anchor_table <- do.call(rbind, lapply(results, `[[`, "anchors"))

print(summary_table)
print(anchor_table)

#write.csv(summary_table, "/Users/nickh/four/DML3_microsynteny_summary.csv", row.names = FALSE)
#write.csv(anchor_table, "/Users/nickh/four/DML3_microsynteny_anchor_pairs.csv", row.names = FALSE)
