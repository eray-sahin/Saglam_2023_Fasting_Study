rm(list = ls())

library(ggplot2)
library(tidyverse)

# Data import and arrangements --------------------------------------------

# metadata
meta <- read.table("Data/ds_metadata.txt", header = T, sep = "\t", stringsAsFactors = T)

# import raw relative abundance table
genus.raw <- read.table("Data/ds_genus_otu_table_final.tsv", sep = "\t", header = T, row.names = 1)

# import OTU-taxa table
genus.otu.taxa <- read.table("Data/ds_genus_taxa_table.tsv", sep = "\t", header = T, row.names = 1)

# replace OTU IDs in raw table with taxonomy names
rownames(genus.raw) <- genus.otu.taxa$Genus[match(rownames(genus.raw), rownames(genus.otu.taxa))]

# transpose raw genus df (re-normalize to get rid of very small changes and make the totals into 100)
genus <- data.frame(t(vegan::decostand(genus.raw, method = "total", MARGIN = 2)*100))

# add time information
genus$meta <- meta$Time[match(rownames(genus), meta$ID)]
genus$meta <- factor(genus$meta, levels = c("before", "after"))

# Preparing top 5 genera table -----------------------------------------------

genus.top5 <- data.frame(matrix(ncol = 5, nrow = nrow(genus)*6))
colnames(genus.top5) <- c("id", "meta", "pair", "genus", "rel.ab")

for (i in 1:nrow(genus)) {
  # id
  genus.top5[((6*i)-5):(6*i), 1] <- rep(rownames(genus)[i], 6)
  # group
  genus.top5[((6*i)-5):(6*i), 2] <- as.character(rep(meta$Time[match(rownames(genus)[i], meta$ID)], 6))
  # pair
  genus.top5[((6*i)-5):(6*i), 3] <- rep(meta$Pairing[match(rownames(genus)[i], meta$ID)], 6)
  # genus
  ordered <- sort(unlist(genus[i, -1833]), decreasing = T)
  genus.top5[((6*i)-5):((6*i)-1), 4] <- names(ordered)[1:5]
  genus.top5[(6*i), 4] <- "Others"
  # relative abundance
  genus.top5[((6*i)-5):((6*i)-1), 5] <- ordered[1:5]
  genus.top5[(6*i), 5] <- sum(ordered[6:length(ordered)])
} 


# png("Genus/Genus_top5.png", units = "in", width = 15, height = 5, res = 800)

ggplot(genus.top5, aes(fill = reorder(genus, rel.ab, FUN = sum), 
                        y = rel.ab, 
                        x = factor(meta, levels = c("before", "after")))) +
  geom_col() +
  facet_grid(.~pair, scale = "free_x") + 
  ggthemes::scale_fill_tableau(palette = "Tableau 20", direction = -1) +
  labs(y = "Relative Abundance (%)", fill = 'Genus') +
  guides(fill = guide_legend(ncol = 2)) +
  theme_light() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 15),
        plot.title = element_blank(),
        legend.position='right',
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"),
        strip.background = element_rect(fill = NA, color = "black", linewidth = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 

# dev.off()
