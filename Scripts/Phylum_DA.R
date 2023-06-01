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

# manual filtering
survived.phyla <- which(colSums(phylum[, -37] < 10^-1) < nrow(phylum)*0.2)
phylum.filtered <- phylum %>% select(all_of(c(names(survived.phyla), "meta")))


# adding B/F ratio
phylum.filtered <- phylum.filtered %>%
  mutate(bf.ratio = Bacteroidetes/Firmicutes) %>%
  relocate(bf.ratio, .before = meta)

# write.table(phylum.filtered,
#             "Phylum/phylum.filtered.ra.txt",
#             sep = "\t", row.names = T, col.names = NA)


# Significance test -------------------------------------------------------

phylum.paired.wilcoxon_p <- vector()
for (i in 1:(ncol(phylum.filtered)-1)) {
  p <- wilcox.test(phylum.filtered[,i] ~ meta, paired = T, exact = F, data = phylum.filtered)$p.value
  phylum.paired.wilcoxon_p[i] <- p
}

colnames(phylum.filtered)[which(phylum.paired.wilcoxon_p < 0.05)]


# creating statistics results for all the filtered phyla ------------------


phylum_filt_all_summary <- data.frame(matrix(ncol = 4, nrow = ncol(phylum.filtered %>% select(!meta))))
colnames(phylum_filt_all_summary) <- c("Before_median",
                                  "After_median",
                                  "After/Before_Log2_Fold_Change", "p.value")
rownames(phylum_filt_all_summary) <- colnames(phylum.filtered %>% select(!meta))

for (i in 1:nrow(phylum_filt_all_summary)) {
  phylum_filt_all_summary$Before_median[i] <- unname(summary(phylum.filtered[phylum.filtered$meta %in% "before", i])[3])
  phylum_filt_all_summary$After_median[i] <- unname(summary(phylum.filtered[phylum.filtered$meta %in% "after", i])[3])
  phylum_filt_all_summary$`After/Before_Log2_Fold_Change`[i] <- log2(phylum_filt_all_summary$After_median[i]/phylum_filt_all_summary$Before_median[i])
  phylum_filt_all_summary$p.value[i] <- phylum.paired.wilcoxon_p[i]
}

# write.table(phylum_filt_all_summary,
#             "Phylum/filtered_phylum_all_stats.tsv",
#             sep = "\t", row.names = T, col.names = NA)


# Plots  ------------------------------------------------------------------


phyla_sig <- phylum.filtered[, c(colnames(phylum.filtered)[which(phylum.paired.wilcoxon_p < 0.05)], "meta")]
phyla_sig$paired <- meta$Pairing[match(rownames(phyla_sig), meta$ID)]

phyla_sig_melted <- reshape2::melt(phyla_sig %>% rownames_to_column(var = "ID"),
                                   id.vars = c("ID", "meta", "paired"))

# adding trendline info to the melted dataframe for plotting

phyla_sig_melted$Trend <- character(48)

for (j in unique(phyla_sig_melted$variable)) {
  # a subset dataframe (df) was created to have a simpler following subtraction
  df <- phyla_sig_melted[which(phyla_sig_melted$variable == j),]
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
    phyla_sig_melted[which(phyla_sig_melted$variable == j), "Trend"][i] <- value
    phyla_sig_melted[which(phyla_sig_melted$variable == j), "Trend"][i+12] <- value
  }
}


# Firmicutes

# png("Phylum/firmicutes.png", units = "in", width = 5, height = 5, res = 800)

ggplot(phyla_sig_melted[which(phyla_sig_melted$variable == "Firmicutes"),], aes_string(x="meta", y="value", fill = "meta")) +
  geom_violin(width=1.0, lwd = 1) +
  geom_boxplot(width=0.15, lwd = 1, color="black", alpha=0.2) +
  geom_point() +
  labs(title = "Firmicutes", x = "Time", y="Relative Abundance (%)") +
  scale_x_discrete(labels=c("before" = "Before", "after" = "After")) +
  scale_fill_manual(name = "Time", values = alpha(c("#EFCA93", "#AFDEA0"), 0.4)) + 
  geom_line(aes(group = paired, color=Trend), size=1, alpha=0.8) +
  scale_color_manual(values = alpha(c("#6666FF","#FF6666"), 0.2)) +
  ggsignif::geom_signif(comparisons = list(c("before", "after")), map_signif_level=function(p) print("p = 0.0038"), size = 1, textsize = 5) +
  theme_light() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        plot.title = element_text(size = 20, face = "italic"),
        legend.position='none') 

# dev.off()

# Proteobacteria

# png("Phylum/proteobacteria.png", units = "in", width = 5, height = 5, res = 800)

ggplot(phyla_sig_melted[which(phyla_sig_melted$variable == "Proteobacteria"),], aes_string(x="meta", y="value", fill = "meta")) +
  geom_violin(width=1.0, lwd = 1) +
  geom_boxplot(width=0.15, lwd = 1, color="black", alpha=0.2) +
  geom_point() +
  labs(title = "Proteobacteria", x = "Time", y="Relative Abundance (%)") +
  scale_x_discrete(labels=c("before" = "Before", "after" = "After")) +
  scale_fill_manual(name = "Time", values = alpha(c("#EFCA93", "#AFDEA0"), 0.4)) + 
  geom_line(aes(group = paired, color=Trend), size=1, alpha=0.8) +
  scale_color_manual(values = alpha(c("#FF6666"), 0.2)) +
  ggsignif::geom_signif(comparisons = list(c("before", "after")), map_signif_level=function(p) print("p = 0.0038"), size = 1, textsize = 5) +
  theme_light() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        plot.title = element_text(size = 20, face = "italic"),
        legend.position='none') 

# dev.off()


