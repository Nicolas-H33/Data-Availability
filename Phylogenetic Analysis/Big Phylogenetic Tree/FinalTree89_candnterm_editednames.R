library(ggtree)
library(ape)
library(dplyr)
library(ggplot2)

setwd("/Users/nickh245/Library/CloudStorage/OneDrive-IowaStateUniversity/Objective 1/Phylogenetic analysis of Gm demeth/PhyloGenetic Trees/RAxML-NG/Final seqs PhyloTrees/Finalseqs(89)/Gm89-candnterm-JTT")
tree <- read.tree("Gm89-candnterm-JTT6.raxml.support")

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
                "G.soja.10G169200.1" = "GsDME1-2", 
                
                "Greenbean.007g223600.1" = "PvROS2",
                "Greenbean.001g186500.1" = "PvROS1",
                "Greenbean.003g032300.1" = "PvDML3",
                
                "Chickpea20137.1" = "CaROS2", 
                "Chickpea17547.1" = "CaROS1",
                "Chickpea08414.1" = "CaDME1",
                "Chickpea07031.1" = "CaDML3",
                
                "Peanut.GS0SHW.1" = "AhROS1-1",
                "Peanut.66FLD0.1" = "AhDME1-1",
                "Peanut.N246S2.1" = "AhROS1-2",
                "Peanut.K8W5NT.1" = "AhROS2-2",
                "Peanut.AM3VD0.1" = "AhROS2-1",
                "Peanut.CUT3DD.2" = "AhDME1-2", 
                "Peanut.R549UJ.2" = "AhDML3-1",
                "Peanut.YBWK9U.1" = "AhDML3-2",
                
                "Redclover.mRNA5898" = "TpROS2",
                "Redclover.mRNA2573" = "TpDME1",
                "Redclover.mRNA33454" = "TpDML3-2",
                "Redclover.mRNA10598" = "TpROS1",
                "Redclover.mRNA35787" = "TpDML3-1",
                
                "Tomato.10T002647.2" = "SlROS1-1",
                "Tomato.03T003535.1" = "SlDML3-1",
                "Tomato.11T000248.1" = "SlDME1",
                "Tomato.09T000303.1" = "SlROS1-2",
                
                "Grape.208s0007g03920.2" = "VvROS1",
                "Grape.213s0074g00450.4" = "VvDME1",
                "Grape.206s0061g01270.1" = "VvDML3",
                
                "Rice.01g11900.1" = "OsDML5-1",
                "Rice.05g37410.1" = "OsDML6-1",
                "Rice.05g37350.1" = "OsDML6-2",
                "Rice.02g29380.1" = "OsDML3-1",
                "Rice.02g29230.1" = "OsDML5-2",
                "Rice.04g28860.1" = "OsDML3-2",
                
                "Apple.MD14G1015300.1" = "MdROS1-1", 
                "Apple.MD12G1017600.1" = "MdROS1-2", 
                "Apple.MD09G1259100" = "MdDME1", 
                "Apple.MD11G1157400" = "MdDML3", 
                
                "Peach.7G118000.1.p" = "PpROS1", 
                "Peach.7G005000.1.p" = "PpDME1-1",
                "Peach.2G089400.1.p" = "PpDME1-2",
                "Peach.6G119100.1.p" = "PpDML3",
                
                "B.cottonwood.006G116000.4" = "PtROS1",
                "B.cottonwood.008G025900.2" = "PtDME1-1",
                "B.cottonwood.010G234400.4" = "PtDME1-2",
                "B.cottonwood.001G150001.7" = "PtDML3",
                
                "Sorghum.004G149800.2" = "SbDML5-2",
                "Sorghum.008G085300.1" = "SbDML5-1",
                "Sorghum.009G155900.1" = "SbDML6",
                "Sorghum.006G224100.1" = "SbDML3",
                
                "C.violacea.0055s0131.1" = "CvROS1",
                "C.violacea.0011s1017.1" = "CvDME1",
                "C.violacea.0011s0895.1" = "CvDML2",
                "C.violacea.0082s0067.1" = "CvDML3",
                
                "B.rapa.E00835.1" = "BrROS1",
                "B.rapa.J02727.1" = "BrDME1",
                "B.rapa.A00413.1" = "BrDML3",
                
                "Corn.eb241240.P002" = "ZmDML5-1",
                "Corn.eb202980.P001" = "ZmDML5-2",
                "Corn.eb289030.P003" = "ZmDML6",
                "Corn.eb241310.P001" = "ZmDML3",
                
                "Cotton.A04G066400.1" = "GhROS1-1(A)",
                "Cotton.A09G136700.1" = "GhROS1-2(A)",
                "Cotton.A01G129700.2" = "GhDME1(A)",
                "Cotton.A05G083800.2" = "GhDML3(A)",
                "Cotton.D04G104500.1" = "GhROS1-1(D)",
                "Cotton.D09G137800.1" = "GhROS1-2(D)",
                "Cotton.D01G119000.1" = "GhDME1(D)",
                "Cotton.D05G086701.1" = "GhDML3(D)",
                
                "Potato.10G024770.1" = "StROS1-1",
                "Potato.09G004240.2" = "StROS1-2",
                "Potato.11G005260.1" = "StDME1",
                "Potato.04G021680.1" = "StDML4-1",
                "Potato.04G017430.1" = "StDML4-2",
                "Potato.03G037240.1" = "StDML3-1")

rooted_tree$tip.label <- ifelse(
  rooted_tree$tip.label %in% names(rename),
  rename[rooted_tree$tip.label],
  rooted_tree$tip.label)

rooted_tree$node.label[rooted_tree$node.label == "Root"] <- ""

#Creation of a data frame for the node labels (numeric values)
node_labels <- data.frame(node = 1:length(rooted_tree$node.label), 
                          label = rooted_tree$node.label,
                          stringsAsFactors = FALSE) 
node_labels$numeric_label <- as.numeric(node_labels$label)

#ggtree plot
p <- ggtree(rooted_tree, layout = "circular", branch.length = "none",  aes(color = as.numeric(label))) + 
  xlim(-1.5, NA) +
  geom_tree(size = 2.6) +
  geom_tiplab(size = 4.5, offset = 0.15, fontface = "bold") + 
  #geom_nodelab(aes(label = label), vjust = 0.8, hjust = 0.02, size = 1, na.rm = TRUE) +
  #geom_treescale(x = 1.8, y = 1, linesize = 0.6, width = 0.5, fontsize = 3.5, color = "black") +
  geom_rootedge(rootedge = 1.5, size = 2.2) +
  geom_cladelab(node = 169, label = "Monocot clades", align = TRUE, offset = 6.05, offset.text = 1, vjust = 3, hjust = -0.15, textcolor = "black", geom = "label", fill = "pink3", barcolor = "pink3", barsize = 2) +
  geom_cladelab(node = 174, label = "", align = TRUE, offset = 6.05, geom = "bar", barcolor = "pink3", barsize = 2) +
  geom_cladelab(node = 93, label = "DML3", align = TRUE, offset = 6.05, offset.text = 1, textcolor = "black", geom = "label", fill = "gold3", barcolor = "gold3", barsize = 2) +
  geom_cladelab(node = 149, label = "DME", align = TRUE, offset = 6.05, offset.text = 2.5, textcolor = "black", geom = "label", fill = "green4", barcolor = "green4", barsize = 2) +
  geom_cladelab(node = 116, label = "ROS/DML2", align = TRUE, offset = 6.05, offset.text = 1, textcolor = "black", geom = "label", fill = "purple2", barcolor = "purple2", barsize = 2) +
  scale_color_gradient(name = "Bootstrap Gradient", low = "blue", high = "red", na.value = "black", guide = guide_colorbar(
    barwidth = 2.5, barheight = 8.5, direction = "vertical", title.position = "top",title.hjust = 0.6)) +
  theme_tree() +
  theme(legend.position = c(0.88, 0.22))
p <- rotate_tree(p, angle = -6)
p <- p + theme(plot.margin = margin(-100, -100, -100, -100))

print(p)
