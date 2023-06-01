rm(list = ls())

library(ggplot2)
library(tidyverse)
library(corrplot)

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

# manual filtering
survived.genera <- which(colSums(genus[, -1833] < 5*10^-2) < nrow(genus)*0.2)
genus.filtered <- genus %>% select(all_of(names(survived.genera)))

# combining metadata and filtered genus abundance

genus.filtered.meta <- genus.filtered %>%
  rownames_to_column(var = "ID") %>%
  left_join(., meta, by = "ID") %>%
  column_to_rownames(var = "ID")

# Removing genera changing with age to exclude from downstream corr --------
age.corr.genera.p <- vector()
last_genus_index_before_age <- grep('Time', colnames(genus.filtered.meta))-1
for (i in 1:last_genus_index_before_age) {
  age.corr.genera.p[i] <- cor.test(genus.filtered.meta$Age, genus.filtered.meta[, i], 
                                   method = "spearman", exact = F)$p.value
}

colnames(genus.filtered.meta)[which(age.corr.genera.p < 0.05)]
age.corr <- which(age.corr.genera.p < 0.05)
genus.filtered.meta.no.age <- genus.filtered.meta[, -age.corr]



# Correlation between genera and nutrients --------------------------------

genus.cor <- genus.filtered.meta.no.age[ , 1:33]
meta.cor <- genus.filtered.meta.no.age[ , 46:70]

genus.feed.cor <- cor(genus.cor, meta.cor, method = "spearman")

genus.feed.cor.p <- matrix(nrow = ncol(genus.cor), ncol = ncol(meta.cor))
rownames(genus.feed.cor.p) <- colnames(genus.cor)
colnames(genus.feed.cor.p) <- colnames(meta.cor)

for (i in 1:ncol(genus.cor)) {
  for (j in 1:ncol(meta.cor)) {
    genus.feed.cor.p[i, j] <- cor.test(genus.cor[, i], meta.cor[, j], method = "spearman", exact = F)$p.value
  }
}


# Correlation Plot --------------------------------------------------------

# illustrating the ones with rho >= 0.5

high.corr <- which(abs(genus.feed.cor) >= 0.5, arr.ind = T)
genus.feed.cor.p[high.corr]
colnames(genus.feed.cor)[high.corr[,2]]


genus.feed.cor_high.corr <- genus.feed.cor[unique(high.corr[,1]), unique(high.corr[,2])]
genus.feed.cor.p_high.corr <- genus.feed.cor.p[unique(high.corr[,1]), unique(high.corr[,2])]

high2 <- which(abs(genus.feed.cor_high.corr) >= 0.5)
genus.feed.cor.p_high.corr[high2]


# png(file="Correlation/genus.feed.cor_high.corr.png", type = "cairo", units = "in", width = 15, height = 20, res = 800)
corrplot(genus.feed.cor_high.corr, 
         p.mat = genus.feed.cor.p_high.corr,
         sig.level = 0.013, # to prevent appearance of ones with abs(rho) < 0.5
         cl.ratio = 0.2, 
         insig = "blank", 
         tl.col = "black",
         tl.srt = 90,
         tl.cex = 3, 
         cl.cex = 3, 
         col = viridis::viridis(100))
# dev.off()
