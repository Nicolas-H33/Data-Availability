library(gggenomes)
library(seqinr)
library(broman)

setwd("/Users/nickh245/Library/CloudStorage/OneDrive-IowaStateUniversity/Objective 1/GmROS & GmDME seqs & guides/Soybean fasta and alignment files/All 5 Gm Whole Genomic Sequence")
WSFile <- read_seqs("All 5 Gm WS.fasta")

#UTR
utr_annotation <- tibble::tibble(
  seq_id = c("GmROS2-1", "GmROS2-1", "GmROS2-1", "GmROS2-1", "GmROS2-1",
             
             "GmROS2-2", "GmROS2-2", "GmROS2-2",
             
             "GmROS1", "GmROS1",
             
             "GmDME1-1", "GmDME1-1",
             
             "GmDME1-2", "GmDME1-2"),
  
  type = c("utr1", "utr2", "utr3", "utr4", "utr5",
           
           "utr1*", "utr2*", "utr3*",
           
           "utr1-", "utr2-",
           
           "utr1#", "utr2#",
           
           "utr1^", "utr2^"),
  
  start = c(1, 250, 421, 538, 16830,
            
            1, 581, 11536,
            
            72, 11395,
            
            191, 11596,
            
            1, 11591),
  
  end = c(164, 324, 454, 707, 17176,
          
          464, 757, 11962,
          
          779, 11765,
          
          1873, 11936,
          
          865, 11927))

#Exon
exon_annotation <- tibble::tibble(
  seq_id = c(rep("GmROS2-1", 19),
             rep("GmROS2-2", 19),
             rep("GmROS1", 20),
             rep("GmDME1-1", 19),
             rep("GmDME1-2", 20)),
  
  type = c(
    #GmROS2-1
    "exon1", "exon2", "exon3", "exon4", "exon5", "exon6", "exon7", "exon8", "exon9", "exon10",
    "exon11", "exon12", "exon13", "exon14", "exon15", "exon16", "exon17", "exon18", "exon19",
    #GmROS2-2
    "exon1*", "exon2*", "exon3*", "exon4*", "exon5*", "exon6*", "exon7*", "exon8*", "exon9*", "exon10*",
    "exon11*", "exon12*", "exon13*", "exon14*", "exon15*", "exon16*", "exon17*", "exon18*", "exon19*",
    #GmROS1
    "exon1-", "exon2-", "exon3-", "exon4-", "exon5-", "exon6-", "exon7-", "exon8-", "exon9-", "exon10-",
    "exon11-", "exon12-", "exon13-", "exon14-", "exon15-", "exon16-", "exon17-", "exon18-", "exon19-", "exon20-",
    #GmDME1-1
    "exon1#", "exon2#", "exon3#", "exon4#", "exon5#", "exon6#", "exon7#", "exon8#", "exon9#", "exon10#",
    "exon11#", "exon12#", "exon13#", "exon14#", "exon15#", "exon16#", "exon17#", "exon18#", "exon19#",
    #GmDME1-2
    "exon1^", "exon2^", "exon3^", "exon4^", "exon5^", "exon6^", "exon7^", "exon8^", "exon9^", "exon10^",
    "exon11^", "exon12^", "exon13^", "exon14^", "exon15^", "exon16^", "exon17^", "exon18^", "exon19^", "exon20^"
  ),
  
  start = c(
    708,1651,4154,5265,5437,6715,6979,7134,7299,7481,7600,7985,8717,8835,9048,9336,9496,9697,16634,
    758,1692,3953,5088,5260,6687,6969,7117,7279,7481,7596,7981,8711,8828,9049,9337,9497,9699,11319,
    780,2057,3591,4545,4720,5725,6114,6386,6541,6717,7088,7217,7672,8474,8590,8848,9138,9303,9492,11184,
    1874,2951,4433,5570,5770,7280,7508,7674,7842,7997,8159,9183,9785,9906,10173,10413,10588,10799,11490,
    866,2092,3181,4831,5920,6134,7642,8112,8291,8459,8614,8779,9125,9786,9910,10181,10417,10593,10802,11487
  ),
  
  end = c(
    1115,4054,4532,5354,6626,6796,7062,7226,7367,7515,7682,8460,8752,8891,9191,9400,9610,9735,16829,
    1165,3864,4331,5177,6596,6768,7052,7209,7347,7515,7678,8456,8746,8884,9192,9401,9611,9737,11535,
    1268,3483,3995,4634,5655,6022,6196,6469,6633,6785,7122,7299,8144,8509,8645,8989,9202,9417,9530,11394,
    2623,4357,4791,5660,7200,7359,7592,7765,7910,8031,8241,9622,9819,9963,10313,10477,10702,10837,11595,
    888,2795,4749,5175,6010,7562,7721,8196,8382,8527,8650,8861,9562,9822,9967,10317,10482,10707,10840,11590
  ))

#Glycosylase domain
glyco_domain <- tibble::tibble(
  seq_id = c(rep("GmROS2-1", 10),
             rep("GmROS2-2", 10),
             rep("GmROS1", 10),
             rep("GmDME1-1", 10),
             rep("GmDME1-2", 10)),
  
  type = c(
    "glyco1","glyco2","glyco3","glyco4","glyco5","glyco6","glyco7","glyco8","glyco9","glyco10",
    "glyco1*","glyco2*","glyco3*","glyco4*","glyco5*","glyco6*","glyco7*","glyco8*","glyco9*","glyco10*",
    "glyco1-","glyco2-","glyco3-","glyco4-","glyco5-","glyco6-","glyco7-","glyco8-","glyco9-","glyco10-",
    "glyco1#","glyco2#","glyco3#","glyco4#","glyco5#","glyco6#","glyco7#","glyco8#","glyco9#","glyco10#",
    "glyco1^","glyco2^","glyco3^","glyco4^","glyco5^","glyco6^","glyco7^","glyco8^","glyco9^","glyco10^"
  ),
  
  start = c(
    4477,5265,5437,6533,6715,6979,7134,7299,7481,7600,
    4276,5088,5260,6503,6687,6969,7117,7279,7481,7596,
    3934,4545,4720,5930,6114,6386,6541,6717,7088,7217,
    4736,5570,5770,7105,7280,7508,7674,7842,7997,8159,
    5118,5920,6134,7467,7642,8112,8291,8459,8614,8779
  ),
  
  end = c(
    4532,5354,5470,6626,6796,7062,7226,7367,7515,7682,
    4331,5177,5293,6596,6768,7052,7209,7347,7515,7678,
    3995,4634,4753,6022,6196,6469,6633,6785,7122,7299,
    4791,5660,5802,7200,7359,7592,7765,7910,8031,8241,
    5175,6010,6164,7562,7721,8196,8382,8527,8650,8859
  ))

#gggenomes plot
Plot2 <- gggenomes(seqs = WSFile, feats = list(exon_annotation, utr_annotation, glyco_domain))
Plot2 +
  geom_seq(linewidth = 1) +
  geom_seq_label(vjust = 1.6, size = 4) +
  geom_feat(aes(color = "Exons"), linewidth = 7.25, data = feats(exon_annotation)) +
  geom_feat(aes(color = "UTR"), linewidth = 3.75, data = feats(utr_annotation)) +
  geom_feat(aes(color = "Glycosylase Domain"), linewidth = 5.25, data = feats(glyco_domain)) +
  scale_color_manual(
    name = "Legend",
    values = c("Exons" = "yellow3", "UTR" = "darkcyan", "Glycosylase Domain" = "firebrick1")
  ) +
  scale_x_continuous(
    limits = c(0, 18500),
    name = "Kilobases",
    breaks = seq(0, 18500, by = 2000),
    labels = scales::label_number(scale = 1e-3),
    expand = c(0, 0)
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 4))) +
  ggtitle("Exon/Intron Structure") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = c(0.875, 0.0825),
    legend.justification = c(0.6, 0),
    legend.box.just = "center",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12)
  )
