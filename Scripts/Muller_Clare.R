# ************************************
# Purpose: Plotting annotated heatmaps
# Date: August 2021
# Author: Salomé Carcy
# ************************************




# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries ####
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(pals)

dir.create("./Heatmaps")

# names(phyloseqobjects) # sanity check

# ***********************
# 2. PREPROCESS DATA ####
# ***********************

## 2.2. Separate fecal & sigmoid samples ####
pseq.root <- subset_samples(
	pseq,
	compartment == 'rhizosphere'
) # 20 samples
pseq.root <- prune_taxa(
	taxa_sums(pseq.root)>0,
	pseq.root
) # remove ASVs that are not present anymore

pseq.bulk <- subset_samples(
	pseq,
	compartment == 'bulk_soil'
) # 20 samples
pseq.bulk <- prune_taxa(
	taxa_sums(pseq.bulk)>0,
	pseq.bulk
) # remove ASVs that are not present anymore

## 2.3. Covariates for heatmap labels ####
color.df <- data.frame(
	compartment = sample_data(pseq)[,'compartment'],
	species = sample_data(pseq)[,'species']
)

color.df <- color.df %>%
	arrange(species, compartment)

# Order of samples
sample.order <- rownames(color.df)
# table(color.df$species) # sanity check

# Colors for heatmap
annotationCol <- list(
	species = c(
		'M. conspicua' = "#f66c42",
		'M. haageana' = "#2995bc",
		'M. lanigera' = "#43b284",
		'M. meissneri' = "#fab255"
	), compartment = c(
		bulk_soil = 'sienna',
		rhizosheath = 'wheat')
)

# *****************************
# 3. HEATMAP PHYLA AS ROWS ####
# *****************************

# # Agglomerate to Phylum level, keeping only ASVs present in at least 2 samples
agglomerate_ps <- function(ps, rank){ #Function to agglomerate by desired taxonomic rank
	res.agg <- ps %>%
		tax_glom(taxrank = rank) %>% # agglomerate at phylum level
		transform_sample_counts(function(x) {x/sum(x)}) %>% # Transform to rel.
		psmelt() # Melt to long format
	return(res.agg)
}

# # Identify phyla present in at least 3 datasets
detect_taxa_ps <- function(ps, rank, fct, cut_off){
	res.agg <- agglomerate_ps(ps, rank)
	res.ls <- res.agg %>%
	# is each phylum present in each dataset (T/F)?
		group_by(.data[[rank]], .data[[fct]]) %>%
		summarize(rank_present = sum(Abundance) > 0, .groups = "drop") %>%
		ungroup %>%
	# in how many datasets is each phylum present (n)?
		group_by(.data[[rank]]) %>%
		count(rank_present) %>%
		filter(rank_present == TRUE) %>%
		filter(n > cut_off) %>%
		ungroup()
	res.ls <- res.ls[[rank]]
	return(res.ls)
}

list.phylum <- detect_taxa_ps(pseq, "Domain", "species", 3)

# # Agglomerate again at Phylum level, but keeping only phyla present in at least 3 datasets

subset_detected <- function(ps, rank, fct, cut_off){
	res.ls <- detect_taxa_ps(ps, rank, fct, cut_off)
	res.agg <- prune_taxa(
		ps@tax_table[, rank] %in% res.ls,
		ps
	)
	res.agg <- res.agg %>%
		tax_glom(taxrank = rank) %>%
		transform_sample_counts(function(x) {x/sum(x)}) %>%
		psmelt()
	return(res.agg)
}

phylum.agg <- subset_detected(pseq, "Phylum", "species", 3)

Get dataframe phylum x samples

phylum.table <- acast(
	phylum.agg %>% filter(Phylum %in% list.phylum),
	Phylum ~ Sample,
	fun.aggregate = sum,
	value.var = 'Abundance'
)

# # Sanity checks
# dim(phylumTable)
# table(is.na(phylumTable))
# table(colSums(phylumTable))
# table(rownames(phylumTable))
# table(rowSums(phylumTable) == 0)
 
# # For coloring, add "pseudocounts"
min(phylum.table[phylum.table>0]) # min is 3.617356e-05
phylum.table[phylum.table == 0] <- 10e-6
 
# # Reorder samples
phylum.table <- phylum.table[, sample.order] # reorder samples

jpeg(
	"./Heatmaps/phylum_heatmap.jpeg",
	height = 3000,
	width = 4000,
	res = 300
)

jpeg(
	"./Heatmaps/corrected_genus_heatmap.jpeg",
	height = 7500,
	width = 4000,
	res = 300
)

jpeg(
	"./Heatmaps/genus_v2_heatmap.jpeg",
	height = 7500,
	width = 4000,
	res = 300
)

genus.table.reorder <- genus.table[clust.order.genus,]

jpeg("./Heatmaps/genus_sampleclust_RF4_delta.jpeg", height = 7500, width = 4000, res = 300)
ph.gn.den <- pheatmap(
	log10(genus.table.reorder),
	color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
	show_rownames = T,
	show_colnames = T,
	angle_col = 90,
	fontsize_row = 15,
	fontsize = 17,
	cluster_rows = F,
	cluster_cols = T,
	#cutree_rows = 2,
	cutree_cols = 2,
#	clustering_method = 'none',
	annotation_col = color.df,
	annotation_colors = annotationCol,
	main = "Hierarchical clustering (Genera as rows)",
	annotation_legend = F
)

grid.text(
	"log10 (rel. freq.)",
	x = 0.92,
	y = 0.94,
	gp = gpar(fontsize = 20, fontface = "bold"),
	rot = 90
)
dev.off()

#Getting the clustering order
clust.order.genus <- ph.gn.den$gtable$grobs[[4]]$label
# **************************************************
# 5. HEATMAP OF DELTAS ####
# **************************************************

ind <- unique(substr(colnames(phylum.table), 1, 9))

new.delta.genus <- matrix(
	nrow = nrow(genus.table),
	ncol = length(ind),
	dimnames = list(
		rownames(genus.table),
		ind
	)
)

for(i in ind){
	new.delta.genus[,i] <- genus.table[,paste0(i, "3")]/genus.table[,paste0(i, "1")]
}

colnames(new.delta.genus) <- substr(ind, 8, 9)

delta.color.df <- data.frame(
	row.names = colnames(tmp.phylum),	
	"species" = rep(levels(color.df$species), each = 5)
)
delta.color.df$species <- as.factor(delta.color.df$species)
delta.annotationCol <- list(species = c(
	'M. conspicua' = "#f66c42",
	'M. haageana' = "#2995bc",
	'M. lanigera' = "#43b284",
	'M. meissneri' = "#fab255"
))

jpeg("./Heatmaps/corrected_genus_delta.jpeg", height = 7500, width = 4000, res = 300)

jpeg("./Heatmaps/genus_FC4_delta.jpeg", height = 7500, width = 4000, res = 300)
new.ph.gn.Ds <- pheatmap(
	scale(log2(new.delta.genus))+(10*asin(delta.genus)), #abs(delta.genus)^(1/3)#
	color = colorRampPalette(brewer.pal(n = 7, name = "BrBG"))(50),
	show_rownames = T,
	show_colnames = T,
	angle_col = 90, # If "show_colnames = T"; this will rotate labels
	fontsize_row = 15,
	fontsize = 20,
	cluster_rows = T,
	cluster_cols = F,
	# cutree_rows = 2,
	# cutree_cols = 2,
	clustering_method = 'ward.D2',
	annotation_col = delta.color.df,
	annotation_colors = delta.annotationCol,
	main = "Pheatmap of Δs (Genera as rows)",
	annotation_legend = F
)
clust.order.genus <- new.ph.gn.Ds$gtable$grobs[[4]]$label

grid.text(
	"Log2(FC)",
	x = 0.92,
	y = 0.93,
	gp = gpar(fontsize = 20, fontface = "bold"),
	rot = 90
)
dev.off()

install.packages("dendextend")
library("dendextend")
row.dend <- as.dendrogram(new.ph.gn.Ds$tree_row)
c.row.dend <- color_branches(
	row.dend,
	k = 4,
	col = c('sienna', 'wheat', 'darkgoldenrod2', 'darkolivegreen2')
)
jpeg("./Heatmaps/row_dendrogram_coloured.jpeg", height = 7500, width = 1000, res = 300)
plot(
	rev(set(c.row.dend, "branches_lwd", 3)),
	horiz = TRUE
)
dev.off()

score.delta.genus <- matrix(
	nrow = nrow(genus.table),
	ncol = length(ind),
	dimnames = list(
		rownames(genus.table),
		ind
	)
)

log2.max <- max(abs(log2(genus.table[,paste0(ind, "3")]/genus.table[,paste0(ind, "1")])))
penalty.max <- max((abs(genus.table[,paste0(ind, "3")]-genus.table[,paste0(ind, "1")])^(1/3)))

for(i in ind){
	RS.i <- genus.table[,paste0(i, "3")]
	BS.i <- genus.table[,paste0(i, "1")]
	log2.i <- log2(RS.i/BS.i)
	penalty.i <- sign(RS.i-BS.i)*(abs(RS.i-BS.i)^(1/3))
	score.delta.genus[,i] <- (log2.i/log2.max)+(penalty.i/penalty.max)
}

jpeg("./Heatmaps/genus_FC6_delta.jpeg", height = 7500, width = 4000, res = 300)
score.ph.gn.Ds <- pheatmap(
	score.delta.genus,
	color = colorRampPalette(brewer.pal(n = 7, name = "BrBG"))(50),
	show_rownames = T,
	show_colnames = F,
	fontsize_row = 15,
	fontsize = 20,
	cluster_rows = T,
	cluster_cols = F,
	# cutree_rows = 2,
	# cutree_cols = 2,
	clustering_method = 'ward.D2',
	annotation_col = delta.color.df,
	annotation_colors = delta.annotationCol,
	main = "Pheatmap of Δs (Genera as rows)",
	annotation_legend = F
)
clust.order.score.genus <- score.ph.gn.Ds$gtable$grobs[[4]]$label

grid.text(
	"Δ-Score",
	x = 0.92,
	y = 0.93,
	gp = gpar(fontsize = 20, fontface = "bold"),
	rot = 90
)
dev.off()

# **************************************************
# 6. HEATMAP OF CORE ASVS ###
# **************************************************

prune_core <- transform_sample_counts(pseq, function(x){x/sum(x)})
prune_core <- prune_taxa(rownames(taxa_core_v2), prune_core)
tmp.naming.vector <- as.vector(prune_core@tax_table[,"Genus"])
tmp.naming.vector[is.na(tmp.naming.vector)] <- as.vector(prune_core@tax_table[is.na(tmp.naming.vector),"Family"])
tmp.naming.vector[is.na(tmp.naming.vector)] <- as.vector(prune_core@tax_table[is.na(tmp.naming.vector),"Order"])
tmp.naming.vector[is.na(tmp.naming.vector)] <- as.vector(prune_core@tax_table[is.na(tmp.naming.vector),"Class"])
tmp.naming.vector[is.na(tmp.naming.vector)] <- as.vector(prune_core@tax_table[is.na(tmp.naming.vector),"Phylum"])
tmp.naming.vector[is.na(tmp.naming.vector)] <- as.vector(prune_core@tax_table[is.na(tmp.naming.vector),"Domain"])
tmp.naming.vector <- paste(tmp.naming.vector, ": ", rownames(prune_core@tax_table))
prune.table <- as.data.frame(prune_core@otu_table)
rownames(prune.table) <- tmp.naming.vector
min(prune.table[prune.table > 0]) # Minimum magnitude is 1e-4
prune.table[prune.table == 0] <- 1e-5
# rm(tmp.naming.vector)

jpeg("./core_pheatmap_v2.jpeg", height = 2500, width = 5000, res = 300)
pheatmap(
	log10(prune.table),
	color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
	show_rownames = T,
	show_colnames = T,
	angle_col = 90,
	fontsize_row = 15,
	fontsize = 20,
	cluster_rows = T,
	cluster_cols = T,
	cutree_rows = 3,
	cutree_cols = 3,
	clustering_method = 'ward.D',
	annotation_col = color.df,
	annotation_colors = annotationCol,
	main = "Pheatmap for the 17 core ASVs",
	annotation_legend = F
)
dev.off()

jpeg("./core_venn_v2.jpeg", height = 2500, width = 3750, res = 300)
core.list
dev.off()

# **************************************************
# 5. HEATMAP OF DELTAS ####
# **************************************************

## I filtered out Archaea by using: pseq <- subset_taxa(pseq, Domain != "Archaea")

library("castor")
library("ape")

phyla.coll <- tax_table(pseq)[,"Phylum"]

# Create a phylum-level tree by finding the MRCA for each phylum
create_phylum_tree <- function(tree, rank){
	unique_taxa <- unique(rank) # Get unique phyla
	taxa_tree <- stree(length(unique_taxa)) # Create a new tree with one tip per phylum
	taxa_tree$tip.label <- unique_taxa
	# Calculate distances between phyla based on their MRCAs
  for (i in 1:(length(unique_phyla)-1)) {
    for (j in (i+1):length(unique_phyla)) {
      # Get ASVs for each phylum
      asvs_i <- names(phyla)[phyla == unique_phyla[i]]
      asvs_j <- names(phyla)[phyla == unique_phyla[j]]
      
      # Find MRCA for each phylum's ASVs
      mrca_i <- get_mrca_of_set(tree, asvs_i)
      mrca_j <- get_mrca_of_set(tree, asvs_j)
      
      # Calculate distance between these MRCAs
      dist <- cophenetic.phylo(tree)[mrca_i, mrca_j]
      
      # Set distance in the new tree
      phylum_tree$edge.length[which(phylum_tree$edge[,2] == j)] <- dist
    }
  }
  
  return(phylum_tree)
}

# Use the function
phylum_tree <- create_phylum_tree(tree, asv_phyla)


# **************************************************
# XXX. DEPRECATED ####
# **************************************************




clust.phylum <- pheatmap(
	tmp.phylum,
	cluster_rows = T,
	cluster_cols = F,
	main = "Δ frequency (rhizosheath − soil) per Phylum"
)

clust_phylum

heat_phylum <- pheatmap( #no clustering hosts
	delta_phylum,
	cluster_rows = T,
	cluster_cols = F,
	main = "Δ frequency (rhizosphere − soil) per phylum",
	color = colorRampPalette(
		c("red", "white", "green")
	)(50),
	breaks = seq(-0.32, 0.32, length.out = 51),
	gaps_col = c(5, 10, 15)
)

heat_phylum

##Averaged deltas
aver_phylum <- long_phylum %>%
	group_by(phylum, site) %>%
	summarise(
		x_delta = mean(delta),
		.groups='drop'
	) %>%
	pivot_wider(
		names_from = site,
		values_from = x_delta,
		values_fill = 0
	) %>% column_to_rownames("phylum")
host_phylum <- pheatmap(
	aver_phylum,
	cluster_rows = T,
	cluster_cols = F,
	main = "Δ frequency (rhizosphere − soil) per phylum"
)

host_phylum