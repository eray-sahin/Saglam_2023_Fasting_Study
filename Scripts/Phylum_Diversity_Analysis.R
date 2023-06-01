rm(list = ls())

library(ggplot2)
library(phyloseq)
library(vegan)
library(tidyverse)

# Data import and arrangements --------------------------------------------

# metadata
meta <- read.table("Data/ds_metadata.txt", header = T, sep = "\t", stringsAsFactors = T) %>%
  column_to_rownames(var = "ID")

# import raw relative abundance table
phylum.raw <- read.table("Data/ds_phylum_otu_table_final.tsv", sep = "\t", header = T, row.names = 1)

# convert relative abundance into counts by multiplying with pseudo number 10^4 to be suitable for phyloseq
phylum.otu <- phylum.raw*10^4

# import OTU-taxa table
phylum.otu.taxa <- read.table("Data/ds_phylum_taxa_table.tsv", sep = "\t", header = T, row.names = 1)


# Create Phyloseq Object --------------------------------------------------

phylum.phyloseq <- phyloseq(
  otu_table(as.matrix(phylum.otu), taxa_are_rows = T),
  sample_data(meta),
  tax_table(as.matrix(phylum.otu.taxa))
)



# Alpha Diversity ---------------------------------------------------------

phylum.alpha.diversity <- estimate_richness(phylum.phyloseq, measures = c("Observed", "Shannon", "Simpson"))

# round into three digits
phylum.alpha.diversity <- round(phylum.alpha.diversity, 3)

# update column names
colnames(phylum.alpha.diversity) <- c("Observed Phyla", "Shannon Index", "Simpson Index")

# add time information 
phylum.alpha.diversity.meta <- phylum.alpha.diversity %>%
  rownames_to_column(var = "ID") %>%
  left_join(., meta %>% 
              rownames_to_column(var = "ID") %>%
              select(ID, Time, Pairing), 
            by = "ID")

phylum.alpha.diversity.meta$Time <- factor(phylum.alpha.diversity.meta$Time, levels = c("before", "after"))

# Paired Wilcoxon ranked sum test -----------------------------------------

phylum_alpha_diversity_p_values <- vector()
for (i in 2:4) {
  p <- wilcox.test(phylum.alpha.diversity.meta[, i] ~ Time, paired = T, exact = F, data = phylum.alpha.diversity.meta)$p.value
  phylum_alpha_diversity_p_values <- append(phylum_alpha_diversity_p_values, p)
}
names(phylum_alpha_diversity_p_values) <- colnames(phylum.alpha.diversity.meta)[2:4]

phylum_alpha_diversity_p_values <- round(phylum_alpha_diversity_p_values, 3)


# Plot --------------------------------------------------------------------

# Melting
phylum.alpha.diversity.melted <- reshape2::melt(phylum.alpha.diversity.meta, id.vars = c("ID", "Pairing", "Time"))
colnames(phylum.alpha.diversity.melted)[4] <- "Alpha_Diversity"
phylum.alpha.diversity.melted$Alpha_Diversity <- factor(phylum.alpha.diversity.melted$Alpha_Diversity, 
                                                        levels = c("Observed Phyla", "Shannon Index", "Simpson Index"))

# adding trendline info to the melted dataframe for plotting

phylum.alpha.diversity.melted$Trend <- character(72)

for (j in unique(phylum.alpha.diversity.melted$Alpha_Diversity)) {
  # a subset dataframe (df) was created to have a simpler following subtraction
  df <- phylum.alpha.diversity.melted[which(phylum.alpha.diversity.melted$Alpha_Diversity == j),]
  for (i in seq(1,12,1)){
    dif <- df$value[i] - df$value[i+12]
    if(dif > 0) {
      value <- "Increasing"
    } else {
      if(dif == 0) {
        value <- "Even"
      } else {
        value <- "Decreasing"
      }}
    phylum.alpha.diversity.melted[which(phylum.alpha.diversity.melted$Alpha_Diversity == j), "Trend"][i] <- value
    phylum.alpha.diversity.melted[which(phylum.alpha.diversity.melted$Alpha_Diversity == j), "Trend"][i+12] <- value
  }
}


phylum.alpha.diversity.melted$Trend <- as.factor(phylum.alpha.diversity.melted$Trend)


# in order to assign a specific colour to each trend category

trend_colors <- c("#6666FF","#727272", "#FF6666")
names(trend_colors) <- c("Decreasing", "Even", "Increasing")

# Plots

for (i in levels(phylum.alpha.diversity.melted$Alpha_Diversity)) {
  # png(paste0("Diversity/Phylum_", i, ".png"), units = "in", width = 5, height = 5, res = 800)
  plot <- ggplot(phylum.alpha.diversity.melted[which(phylum.alpha.diversity.melted$Alpha_Diversity == i),], 
                 aes_string(x = "Time", y = "value", fill = "Time")) +
    geom_violin(width = 1.0, linewidth = 1) +
    geom_boxplot(width = 0.15, linewidth = 1, color = "black", alpha = 0.2) +
    geom_point() +
    labs(title = i, x = "Time", y = "Alpha Diversity") +
    scale_x_discrete(labels = c("before" = "Before", "after" = "After")) +
    scale_fill_manual(name = "Time", values = alpha(c("#EFCA93", "#AFDEA0"), 0.4)) + 
    geom_line(aes(group = Pairing, color = Trend), size=1, alpha=0.8) +
    scale_color_manual(values = trend_colors) +
    ggsignif::geom_signif(comparisons = list(c("before", "after")), 
                          map_signif_level=function(p) print(paste0("p = ", 
                                                                    phylum_alpha_diversity_p_values[names(phylum_alpha_diversity_p_values) == i])), 
                          size = 1, textsize = 5) +
    theme_light() +
    theme(axis.line = element_line(colour = "grey"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          plot.title = element_text(size = 20),
          legend.position = 'none')
  
  print(plot)
  
  # dev.off()
  
}


# Beta Diversity ----------------------------------------------------------

# Bray-Curtis
phylum_bray_dist = phyloseq::distance(phylum.phyloseq, method = "bray")
phylum_ordination = ordinate(phylum.phyloseq, method = "PCoA", distance = phylum_bray_dist)
phylum_axis1.2 <- as.data.frame(phylum_ordination$vectors[,1:2])
phylum_axis1.2$rel_eigen <- phylum_ordination$values$Relative_eig
colnames(phylum_axis1.2) <- c("PC1", "PC2", "rel_eigen")
phylum_axis1.2 <- phylum_axis1.2 %>% rownames_to_column("ID")
phylum_axis1.2$meta <- meta$Time[match(phylum_axis1.2$ID, rownames(meta))]
phylum_axis1.2$meta <- factor(phylum_axis1.2$meta, levels = c("before", "after"))

# PERMANOVA
set.seed(123)
adonis2(phylum_bray_dist ~ meta, permutations = 9999, data = phylum_axis1.2)

meta_colors <- c("#EFCA93", "#AFDEA0")
names(meta_colors) <- levels(phylum_axis1.2$meta)

# png("Diversity/phylum_beta_diversity.png", units = "in", width = 7.5, height = 5, res = 800)

phylum_axis1.2 %>%
  ggplot(aes(x = PC1, y = PC2, colour = meta)) +
  geom_point(size = 3) + 
  stat_ellipse(linewidth = 2) + 
  scale_color_manual(name = "Time", labels = c("Before", "After"), values = meta_colors) +
  labs(x=paste0("PC1 (", round(phylum_axis1.2$rel_eigen[1] * 100, digits = 2), "%)"), 
       y=paste0("PC2 (", round(phylum_axis1.2$rel_eigen[2] * 100, digits = 2), "%)"),
       title = "Phylum - PCoA (Bray-Curtis)") + theme(axis.line = element_line(colour = "grey"),
                                                     panel.grid.major = element_blank(),
                                                     panel.grid.minor = element_blank(),
                                                     panel.border = element_blank(),
                                                     panel.background = element_blank(),
                                                     axis.title = element_text(size = 15),
                                                     axis.text.x = element_text(size = 15),
                                                     axis.text.y = element_text(size = 15),
                                                     plot.title = element_text(size = 20),
                                                     legend.text = element_text(size = 15),
                                                     legend.title = element_text(size = 15),
                                                     strip.text = element_text(size=15),
                                                     legend.position='right') +
  # Adding text box for PERMANOVA result
  annotate("text", label = "F = 6.353, p = 0.0014", x = -0.27, y = -0.055, size = 5)

# dev.off()
