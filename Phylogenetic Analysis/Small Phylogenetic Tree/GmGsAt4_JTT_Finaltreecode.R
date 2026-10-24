library(ggtree)
library(ape)
library(dplyr)

setwd("/Users/nickh245/Library/CloudStorage/OneDrive-IowaStateUniversity/Objective 1/Phylogenetic analysis of Gm demeth/PhyloGenetic Trees/RAxML-NG/G.max, G.soja, and A.thaliana PhyloTrees/GmGsAt4-candnterm-JTT")
tree <- read.tree("GmGsAt4-candnterm-JTT6.raxml.support")

outgroup_gym <- "C.japonica"
rooted_tree <- root(tree, outgroup = outgroup_gym, resolve.root = TRUE)

rename <- c("GmROS1" = "GmROS2-1",
            "GmROS2" = "GmROS2-2",
            "GmROS3" = "GmROS1",
            "GmDME1" = "GmDME1-1",
            "GmDME2" = "GmDME1-2",
            "G.soja.10G060200.1" = "GsROS2-1",
            "G.soja.20G152600.1" = "GsDME1-1",
            "G.soja.13G117300.1" = "GsROS2-2",
            "G.soja.03G152500.1" = "GsROS1",
            "G.soja.10G169200.1" = "GsDME1-2")

rooted_tree$tip.label <- ifelse(
  rooted_tree$tip.label %in% names(rename),
  rename[rooted_tree$tip.label],
  rooted_tree$tip.label)

rooted_tree$node.label[rooted_tree$node.label == "Root"] <- ""

#ggtree plot
p <- ggtree(rooted_tree) + 
  geom_tree(size = 0.7) +
  geom_tiplab(size = 3.4) + 
  geom_nodelab(aes(label = label), vjust = 0.5, hjust = -0.03, size = 3.15, na.rm = TRUE) +
  geom_treescale(x = 2.4, y = 2.4, width = 1, fontsize = 3.5, color = "black") +
  geom_rootedge(rootedge = 0.2) +
  theme_tree() +
  theme(legend.position = c(0.825, 0.37))

print(p)
