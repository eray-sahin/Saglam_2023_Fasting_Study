rm(list = ls())

library(ggplot2)
library(tidyverse)

# Data import and arrangements --------------------------------------------

# metadata
meta <- read.table("Data/ds_metadata.txt", header = T, sep = "\t", stringsAsFactors = T)

# import raw relative abundance table
phylum.raw <- read.table("Data/ds_phylum_otu_table_final.tsv", sep = "\t", header = T, row.names = 1)

# import OTU-taxa table
phylum.otu.taxa <- read.table("Data/ds_phylum_taxa_table.tsv", sep = "\t", header = T, row.names = 1)

# replace OTU IDs in raw table with taxonomy names
rownames(phylum.raw) <- phylum.otu.taxa$Phylum[match(rownames(phylum.raw), rownames(phylum.otu.taxa))]

# # transpose raw phylum df (re-normalize to get rid of very small changes and make the totals into 100)
phylum <- data.frame(t(vegan::decostand(phylum.raw, method = "total", MARGIN = 2)*100))

# add time information
phylum$meta <- meta$Time[match(rownames(phylum), meta$ID)]
phylum$meta <- factor(phylum$meta, levels = c("before", "after"))


# Preparing top 3 phyla table -----------------------------------------------

phylum.top3 <- data.frame(matrix(ncol = 5, nrow = nrow(phylum)*4))
colnames(phylum.top3) <- c("id", "meta", "pair", "phylum", "rel.ab")

for (i in 1:nrow(phylum)) {
  # id
  phylum.top3[((4*i)-3):(4*i), 1] <- rep(rownames(phylum)[i], 4)
  # group
  phylum.top3[((4*i)-3):(4*i), 2] <- as.character(rep(meta$Time[match(rownames(phylum)[i], meta$ID)], 4))
  # pair
  phylum.top3[((4*i)-3):(4*i), 3] <- rep(meta$Pairing[match(rownames(phylum)[i], meta$ID)], 4)
  # phylum
  ordered <- sort(unlist(phylum[i, -37]), decreasing = T)
  phylum.top3[((4*i)-3):((4*i)-1), 4] <- names(ordered)[1:3]
  phylum.top3[(4*i), 4] <- "Others"
  # relative abundance
  phylum.top3[(4*i)-3, 5] <- ordered[1]
  phylum.top3[(4*i)-2, 5] <- ordered[2]
  phylum.top3[(4*i)-1, 5] <- ordered[3]
  phylum.top3[(4*i), 5] <- sum(ordered[4:length(ordered)])
} 



# Plot --------------------------------------------------------------------


# png("Phylum/Phylum_top3.png", units = "in", width = 15, height = 5, res = 800)

ggplot(phylum.top3, aes(fill=reorder(phylum, rel.ab, FUN = sum), 
                       y=rel.ab, 
                       x=factor(meta, levels = c("before", "after")))) + 
  geom_col() +
  facet_grid(.~pair, scale = "free_x") + 
  ggthemes::scale_fill_tableau(palette = "Tableau 20", direction = -1) +
  labs(y="Relative Abundance (%)", fill='Phylum') +
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

