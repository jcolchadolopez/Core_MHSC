setwd("D:/CleanData")

library("phyloseq")
library("vegan")
rseq <- readRDS("phyloseq_rseq.rds")

library("igraph")
#library("devtools")
#install_github("zdk123/SpiecEasi")
library("SpiecEasi")
library("ggplot2")
library("ggraph")

#Filtering low-abundance
MR_rseq <- subset_samples(rseq, compartment == 2)
filt_995 <- prune_taxa(
	taxa_sums(MR_rseq) > quantile(taxa_sums(MR_rseq), probs = 0.995),
	MR_rseq
)

mod_995 <- spiec.easi(filt_995,
		method = "mb",
		lambda.min.ratio = 0.003,
		nlambda = 100
	)

g_995 <- as.matrix(symBeta(getOptBeta(mod_995))) %>%
		graph_from_adjacency_matrix(
			mode = "undirected",
			weighted = T,
			diag = F
		)
plot(g_995)

V(g_995)$label <- as.character(tax_table(filt_995)[V(g_995), "Genus"])
tmp.tax <- tax_table(filt_995)[V(g_995), "Genus"]
tmp.tax <- intersect(rownames(tmp.tax), rownames(taxa_core))
tmp.tax <- tax_table(filt_995)[tmp.tax,]
#taxa_core is from the Venn.R script

tmp.df <- as.data.frame(tax_table(filt_995)[V(g_995),])
tmp.df$Final <- rep("", nrow(tmp.df))
for(i in 1:nrow(tmp.df)){
	tmp.row <- tmp.df[i,]
#	print(tmp.row)
	if(!is.na(tmp.row$Species) && !grepl("^uncultured_", tmp.row$Species)){
		tmp.df$Final[i] <- tmp.row$Species
		next
	} else if(!is.na(tmp.row$Genus)){
		tmp.df$Final[i] <- paste0(tmp.row$Genus, " sp.")
		next
	} else if(!is.na(tmp.row$Familia)){
		tmp.df$Final[i] <- paste0("f_", tmp.row$Familia)
		next	
	} else if(!is.na(tmp.row$Ordo)){
		tmp.df$Final[i] <- paste0("o_", tmp.row$Ordo)
		next
	} else if(!is.na(tmp.row$Classis)){
		tmp.df$Final[i] <- paste0("c_", tmp.row$Classis)
		next
	} else if(!is.na(tmp.row$Phylum)){
		tmp.df$Final[i] <- paste0("p_", tmp.row$Phylum)
		next
	} else {
		tmp.df$Final[i] <- "Unknown bacteria"
	}
}
tmp.df <- tmp.df[-c(1:7)]

tmp.tax <- as.data.frame(tmp.tax)
tmp.tax$Final <- rep("", nrow(tmp.tax))
for(i in 1:nrow(tmp.tax)){
	tmp.row <- tmp.tax[i,]
#	print(tmp.row)
	if(!is.na(tmp.row$Species) && !grepl("^uncultured_", tmp.row$Species)){
		tmp.tax$Final[i] <- tmp.row$Species
		next
	} else if(!is.na(tmp.row$Genus)){
		tmp.tax$Final[i] <- paste0(tmp.row$Genus, " sp.")
		next
	} else if(!is.na(tmp.row$Familia)){
		tmp.tax$Final[i] <- paste0("f_", tmp.row$Familia)
		next	
	} else if(!is.na(tmp.row$Ordo)){
		tmp.tax$Final[i] <- paste0("o_", tmp.row$Ordo)
		next
	} else if(!is.na(tmp.row$Classis)){
		tmp.tax$Final[i] <- paste0("c_", tmp.row$Classis)
		next
	} else if(!is.na(tmp.row$Phylum)){
		tmp.tax$Final[i] <- paste0("p_", tmp.row$Phylum)
		next
	} else {
		tmp.tax$Final[i] <- "Unknown bacteria"
	}
}
tmp.tax <- tmp.tax[-c(1:7)]

tmp.names <- which(rownames(tmp.df) %in% rownames(tmp.tax), arr.ind = T)

set.seed(995)
plot_995 <- ggraph(g_995, layout = "fr")+
		geom_edge_link(aes(width = weight), color = "green")+
		geom_node_point(aes(size = degree(g_995)), color = "red")+
		geom_node_text(aes(
			label = tmp.df$Final,
			fontface = ifelse(
				rownames(tmp.df) %in% rownames(tmp.df)[tmp.names],
				"bold.italic",
				"italic"
			)
		), repel=T, size =4)
plot_995

#Less constrictive
filt_99 <- prune_taxa(
	taxa_sums(MR_rseq) > quantile(taxa_sums(MR_rseq), probs = 0.99),
	MR_rseq
)

mod_99 <- spiec.easi(filt_99,
		method = "mb",
		lambda.min.ratio = 0.1,
		nlambda = 10
	)

g_99 <- as.matrix(symBeta(getOptBeta(mod_99))) %>%
		graph_from_adjacency_matrix(
			mode = "undirected",
			weighted = T,
			diag = F
		)

V(g_99)$label <- as.character(tax_table(filt_99)[V(g_99), "Genus"])
tmp.tax <- tax_table(filt_99)[V(g_99), "Genus"]
tmp.tax <- intersect(rownames(tmp.tax), rownames(taxa_core))
tmp.tax <- tax_table(filt_99)[tmp.tax,]
#taxa_core is from the Venn.R script

tmp.df <- as.data.frame(tax_table(filt_99)[V(g_99),])
tmp.df$Final <- rep("", nrow(tmp.df))
for(i in 1:nrow(tmp.df)){
	tmp.row <- tmp.df[i,]
#	print(tmp.row)
	if(!is.na(tmp.row$Species) && !grepl("^uncultured_", tmp.row$Species)){
		tmp.df$Final[i] <- tmp.row$Species
		next
	} else if(!is.na(tmp.row$Genus)){
		tmp.df$Final[i] <- paste0(tmp.row$Genus, " sp.")
		next
	} else if(!is.na(tmp.row$Familia)){
		tmp.df$Final[i] <- paste0("f_", tmp.row$Familia)
		next	
	} else if(!is.na(tmp.row$Ordo)){
		tmp.df$Final[i] <- paste0("o_", tmp.row$Ordo)
		next
	} else if(!is.na(tmp.row$Classis)){
		tmp.df$Final[i] <- paste0("c_", tmp.row$Classis)
		next
	} else if(!is.na(tmp.row$Phylum)){
		tmp.df$Final[i] <- paste0("p_", tmp.row$Phylum)
		next
	} else {
		tmp.df$Final[i] <- "Unknown bacteria"
	}
}
tmp.df <- tmp.df[-c(1:7)]

tmp.tax <- as.data.frame(tmp.tax)
tmp.tax$Final <- rep("", nrow(tmp.tax))
for(i in 1:nrow(tmp.tax)){
	tmp.row <- tmp.tax[i,]
#	print(tmp.row)
	if(!is.na(tmp.row$Species) && !grepl("^uncultured_", tmp.row$Species)){
		tmp.tax$Final[i] <- tmp.row$Species
		next
	} else if(!is.na(tmp.row$Genus)){
		tmp.tax$Final[i] <- paste0(tmp.row$Genus, " sp.")
		next
	} else if(!is.na(tmp.row$Familia)){
		tmp.tax$Final[i] <- paste0("f_", tmp.row$Familia)
		next	
	} else if(!is.na(tmp.row$Ordo)){
		tmp.tax$Final[i] <- paste0("o_", tmp.row$Ordo)
		next
	} else if(!is.na(tmp.row$Classis)){
		tmp.tax$Final[i] <- paste0("c_", tmp.row$Classis)
		next
	} else if(!is.na(tmp.row$Phylum)){
		tmp.tax$Final[i] <- paste0("p_", tmp.row$Phylum)
		next
	} else {
		tmp.tax$Final[i] <- "Unknown bacteria"
	}
}
tmp.tax <- tmp.tax[-c(1:7)]

tmp.names <- which(rownames(tmp.df) %in% rownames(tmp.tax), arr.ind = T)

set.seed(99)
g_99_pos <- subgraph.edges(g_99, E(g_99)[weight > 0])
g_99_neg <- subgraph.edges(g_99, E(g_99)[weight < 0])
E(g_99_neg)$weight <- abs(E(g_99_neg)$weight)
g_99_abs <- g_99
E(g_99_abs)$weight <- abs(E(g_99_abs)$weight)
set.seed(99)
plot_99 <- ggraph(g_99_abs, layout = "fr")+
		geom_edge_link(aes(width = E(g_99_abs)$weight), color = "green")+
		geom_node_point(aes(size = degree(g_99_abs)), color = "red")+
		geom_node_text(aes(
			label = tmp.df$Final,
			fontface = ifelse(
				rownames(tmp.df) %in% rownames(tmp.df)[tmp.names],
				"bold.italic",
				"italic"
			)
		), repel=T, size =4)+
		ggtitle("Network topology \n (occurrence probability > 0.99)")+
		theme(plot.title = element_text(
			face = "bold",
			size = 17,
			hjust = 0.5
			)
		)
plot_99

set.seed(99)
neg_99 <- ggraph(g_99_neg, layout = "fr")+
		geom_edge_link(aes(width = weight), color = "red")+
		geom_node_point(aes(size = degree(g_99_neg)), color = "blue")+
		geom_node_text(aes(label = name), repel=T, size =4)+
		ggtitle("Network topology \n (occurrence probability > 0.99)")+
		theme(plot.title = element_text(
			face = "bold",
			size = 17,
			hjust = 0.5
			)
		)
neg_99


degree_dist <- degree(g_99)
betweenness_dist <- betweenness(g_99_abs)
par(mfrow = c(2,1))
plot(
	table(degree_dist),
	type = 'b',
	main = "Plot of degrees per node",
	xlab = substitute(paste(bold("Degrees"))),
	ylab = substitute(paste(bold("Frequency"))),
	cex = 1.2, pch = 19, cex.axis = 1, cex.lab = 1.3, cex.main = 1.5)
)
abline(v = median(degree_dist), col = "red", lwd = 3)
plot(
	table(betweenness_dist),
	type = 'b',
	main = "Plot of betweenness per node",
	xlab = substitute(paste(bold("betweeneess distance"))),
	ylab = substitute(paste(bold("frequency"))),
	cex = 1.2, pch = 19, cex.axis = 1, cex.lab = 1.3, cex.main = 1.5)
)
abline(v = quantile(betweenness_dist, probs = 0.75), col = "red", lwd = 3)

clustering_coeff <- transitivity(g_995, type = "local")