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

# manual filtering
survived.genera <- which(colSums(genus[, -1833] < 5*10^-2) < nrow(genus)*0.2)
genus.filtered <- genus %>% select(all_of(c(names(survived.genera), "meta")))

# write.table(genus.filtered,
#             "Genus/genus.filtered.ra.txt",
#             sep = "\t", row.names = T, col.names = NA)

# Significance test -------------------------------------------------------

paired.wilcoxon_p <- vector()
for (i in 1:(ncol(genus.filtered)-1)){
  p <- wilcox.test(genus.filtered[,i] ~ meta, paired = T, exact = F, data = genus.filtered)$p.value
  paired.wilcoxon_p[i] <- p
}

colnames(genus.filtered)[which(paired.wilcoxon_p < 0.05)]

# creating statistics results for all the filtered genera ------------------

genus_filt_all_summary <- data.frame(matrix(ncol = 4, nrow = ncol(genus.filtered %>% select(!meta))))
colnames(genus_filt_all_summary) <- c("Before_median",
                                       "After_median",
                                       "After/Before_Log2_Fold_Change", "p.value")
rownames(genus_filt_all_summary) <- colnames(genus.filtered %>% select(!meta))

for (i in 1:nrow(genus_filt_all_summary)) {
  genus_filt_all_summary$Before_median[i] <- unname(summary(genus.filtered[genus.filtered$meta %in% "before", i])[3])
  genus_filt_all_summary$After_median[i] <- unname(summary(genus.filtered[genus.filtered$meta %in% "after", i])[3])
  genus_filt_all_summary$`After/Before_Log2_Fold_Change`[i] <- log2(genus_filt_all_summary$After_median[i]/genus_filt_all_summary$Before_median[i])
  genus_filt_all_summary$p.value[i] <- paired.wilcoxon_p[i]
}

# write.table(genus_filt_all_summary,
#             "Genus/filtered_genus_all_stats.tsv",
#             sep = "\t", row.names = T, col.names = NA)


# Plots  ------------------------------------------------------------------

genera_sig <- genus.filtered[, c(colnames(genus.filtered)[which(paired.wilcoxon_p < 0.05)], "meta")]
genera_sig$paired <- meta$Pairing[match(rownames(genera_sig), meta$ID)]

genera_sig_melted <- reshape2::melt(genera_sig %>% rownames_to_column(var = "ID"),
               id.vars = c("ID", "meta", "paired"))

# adding trendline info to the melted dataframe for plotting

genera_sig_melted$Trend <- character(nrow(genera_sig_melted))

for (j in unique(genera_sig_melted$variable)) {
  # a subset dataframe (df) was created to have a simpler following subtraction
  df <- genera_sig_melted[which(genera_sig_melted$variable == j),]
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
    genera_sig_melted[which(genera_sig_melted$variable == j), "Trend"][i] <- value
    genera_sig_melted[which(genera_sig_melted$variable == j), "Trend"][i+12] <- value
  }
}

genera_sig_melted$Trend <- as.factor(genera_sig_melted$Trend)

# in order to assing a specific colour to each trend category

trend_colors <- c("#6666FF","#727272", "#FF6666")
names(trend_colors) <- c("Decreasing", "Even", "Increasing")

# Plots of 9 significantly changing genera

for (i in colnames(genera_sig)[1:9]) {
  
  # png(paste0("Genus/", i, "_plot.png"), units = "in", width = 5, height = 5, res = 800)
  p <- ggplot(genera_sig_melted[which(genera_sig_melted$variable == i),], 
                           aes_string(x="meta", y="value", fill = "meta")) +
    geom_violin(width=1.0, lwd = 1) +
    geom_boxplot(width=0.15, lwd = 1, color="black", alpha=0.2) +
    geom_point() +
    labs(title = i, x = "Time", y="Relative Abundance (%)") +
    scale_x_discrete(labels=c("before" = "Before", "after" = "After")) +
    scale_fill_manual(name = "Time", values = alpha(c("#EFCA93", "#AFDEA0"), 0.4)) + 
    geom_line(aes(group = paired, color = Trend), size = 1, alpha = 0.8) +
    scale_color_manual(values = trend_colors) +
    ggsignif::geom_signif(comparisons = list(c("before", "after")), 
                          map_signif_level=function(p) print(paste0("p = ", 
                                                                    round(genus_filt_all_summary$p.value[match(i, rownames(genus_filt_all_summary))], 3))), 
                          size = 1, textsize = 5) +
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
  
  print(p)
  
  # dev.off()
  
}
