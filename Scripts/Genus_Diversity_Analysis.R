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
genus.raw <- read.table("Data/ds_genus_otu_table_final.tsv", sep = "\t", header = T, row.names = 1)

# convert relative abundance into counts by multiplying with pseudo number 10^4 to be suitable for phyloseq
genus.otu <- genus.raw*10^4

# import OTU-taxa table
genus.otu.taxa <- read.table("Data/ds_genus_taxa_table.tsv", sep = "\t", header = T, row.names = 1)


# Create Phyloseq Object --------------------------------------------------

genus.phyloseq <- phyloseq(
  otu_table(as.matrix(genus.otu), taxa_are_rows = T),
  sample_data(meta),
  tax_table(as.matrix(genus.otu.taxa))
)


# Alpha Diversity ---------------------------------------------------------

genus.alpha.diversity <- estimate_richness(genus.phyloseq, measures = c("Observed", "Shannon", "Simpson"))

# round into three digits
genus.alpha.diversity <- round(genus.alpha.diversity, 3)

# update column names
colnames(genus.alpha.diversity) <- c("Observed Genera", "Shannon Index", "Simpson Index")

# add time information 
genus.alpha.diversity.meta <- genus.alpha.diversity %>%
  rownames_to_column(var = "ID") %>%
  left_join(., meta %>% 
              rownames_to_column(var = "ID") %>%
              select(ID, Time, Pairing), 
            by = "ID")

genus.alpha.diversity.meta$Time <- factor(genus.alpha.diversity.meta$Time, levels = c("before", "after"))

# Paired Wilcoxon ranked sum test -----------------------------------------

genus_alpha_diversity_p_values <- vector()
for (i in 2:4) {
  p <- wilcox.test(genus.alpha.diversity.meta[, i] ~ Time, paired = T, exact = F, data = genus.alpha.diversity.meta)$p.value
  genus_alpha_diversity_p_values <- append(genus_alpha_diversity_p_values, p)
}
names(genus_alpha_diversity_p_values) <- colnames(genus.alpha.diversity.meta)[2:4]

genus_alpha_diversity_p_values <- round(genus_alpha_diversity_p_values, 3)


# Plot --------------------------------------------------------------------

# Melting
genus.alpha.diversity.melted <- reshape2::melt(genus.alpha.diversity.meta, id.vars = c("ID", "Pairing", "Time"))
colnames(genus.alpha.diversity.melted)[4] <- "Alpha_Diversity"
genus.alpha.diversity.melted$Alpha_Diversity <- factor(genus.alpha.diversity.melted$Alpha_Diversity, 
                                                        levels = c("Observed Genera", "Shannon Index", "Simpson Index"))

# adding trendline info to the melted dataframe for plotting

genus.alpha.diversity.melted$Trend <- character(72)

for (j in unique(genus.alpha.diversity.melted$Alpha_Diversity)) {
  # a subset dataframe (df) was created to have a simpler following subtraction
  df <- genus.alpha.diversity.melted[which(genus.alpha.diversity.melted$Alpha_Diversity == j),]
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
    genus.alpha.diversity.melted[which(genus.alpha.diversity.melted$Alpha_Diversity == j), "Trend"][i] <- value
    genus.alpha.diversity.melted[which(genus.alpha.diversity.melted$Alpha_Diversity == j), "Trend"][i+12] <- value
  }
}


genus.alpha.diversity.melted$Trend <- as.factor(genus.alpha.diversity.melted$Trend)


# in order to assign a specific colour to each trend category

trend_colors <- c("#6666FF","#727272", "#FF6666")
names(trend_colors) <- c("Decreasing", "Even", "Increasing")

# Plots

for (i in levels(genus.alpha.diversity.melted$Alpha_Diversity)) {
  # png(paste0("Diversity/Genus_", i, ".png"), units = "in", width = 5, height = 5, res = 800)
  plot <- ggplot(genus.alpha.diversity.melted[which(genus.alpha.diversity.melted$Alpha_Diversity == i),], 
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
                                                                    genus_alpha_diversity_p_values[names(genus_alpha_diversity_p_values) == i])), 
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
genus_bray_dist = phyloseq::distance(genus.phyloseq, method = "bray")
genus_ordination = ordinate(genus.phyloseq, method = "PCoA", distance = genus_bray_dist)
genus_axis1.2 <- as.data.frame(genus_ordination$vectors[,1:2])
genus_axis1.2$rel_eigen <- genus_ordination$values$Relative_eig
colnames(genus_axis1.2) <- c("PC1", "PC2", "rel_eigen")
genus_axis1.2 <- genus_axis1.2 %>% rownames_to_column("ID")
genus_axis1.2$meta <- meta$Time[match(genus_axis1.2$ID, rownames(meta))]
genus_axis1.2$meta <- factor(genus_axis1.2$meta, levels = c("before", "after"))

# PERMANOVA
set.seed(123)
adonis2(genus_bray_dist ~ meta, permutations = 9999, data = genus_axis1.2)

meta_colors <- c("#EFCA93", "#AFDEA0")
names(meta_colors) <- levels(genus_axis1.2$meta)

# png("Diversity/genus_beta_diversity.png", units = "in", width = 7.5, height = 5, res = 800)

genus_axis1.2 %>%
  ggplot(aes(x = PC1, y = PC2, colour = meta)) +
  geom_point(size = 3) + 
  stat_ellipse(linewidth = 2) + 
  scale_color_manual(name = "Time", labels = c("Before", "After"), values = meta_colors) +
  labs(x=paste0("PC1 (", round(genus_axis1.2$rel_eigen[1] * 100, digits = 2), "%)"), 
       y=paste0("PC2 (", round(genus_axis1.2$rel_eigen[2] * 100, digits = 2), "%)"),
       title = "Genus - PCoA (Bray-Curtis)") + theme(axis.line = element_line(colour = "grey"),
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
  annotate("text", label = "F = 1.64, p = 0.082", x = -0.3, y = -0.28, size = 5)

# dev.off()
