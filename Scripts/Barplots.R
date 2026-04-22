#### ESTE SCRIPT EXPLICA EL PROCESO PARA REALIZAR LOS BARPLOTS ####
setwd("D:/CleanData")
library("microbiome")
library("microbiomeutilities")
library("phyloseq")
library("tidyverse")
library("vegan")
library("dplyr")

#### REQUISITOS ####
# Un objeto `phyloseq` con mínimo una tabla de taxonomía y una de abundancias

#Files to load:
biom <- read.table(
	"unmocked.tsv",
	header=F,
	sep='\t'
)

meta <- read.csv(
	"METACONST.csv",
	header=T,
	row.names=1
)

#Corrections to adjust data format
rownames(meta)[23] <- "sample-341"
colnames(biom) <- c("ASV", rownames(meta))
rownames(biom) <- biom[,1]

biom <- biom[,-1]
biom <- as.data.frame(
	t(
		as.matrix(biom)
	)
)

meta <- meta[,-c(3:5)]
biom <- biom[order(rownames(biom)),]

meta$species <- as.factor(meta$species)
meta$compartment <- as.factor(meta$compartment)
meta$treatment <- as.factor(paste0(meta$species, "_", meta$compartment))
levels(meta$treatment) <- c("M. conspicua BS", "M. conspicua MR",
					"M. haageana BS", "M. haagean MR",
					"M. lanigera BS", "M. lanigera MR",
					"M. meissneri BS", "M. meissneri MR"
	)
taxa <- unique(biom)

###Validating data function###
validating.names <- function(x, taxrank, NearestName = TRUE, verbose = TRUE){
	pos.col <- which(colnames(tax_table(x))==taxrank)
	tax.vex <- as.vector(tax_table(x)[,pos.col])
	na_lox <- which(is.na(tax.vex))
	if(length(na_lox) == 0){
		if(verbose) print("All taxa have valid names")
		return(tax.vex)
	} else {
		if(verbose) {
			print(
				paste0(
					"Converting 'NA's to 'Unknown' at positions: ",
					paste(na_lox, collapse = ", ")
				)
			)
		}
		for(i in na_lox){
			if(NearestName){
				current_col <- pos.col - 1
				H.rank <- NA
				while(current_col > 1 && is.na(H.rank)){
					H.rank <- tax_table(x)[i, current_col]
					current_col <- current_col - 1
				}
				H.name <- ifelse(is.na(H.rank), "",
					paste0("_", H.rank))
			} else {
				H.name <- ""
			}
			new_name <- paste0("Unknown", H.name)
			tax_table(x)[i, pos.col] <- new_name
			if(verbose) print(paste0(i, " : ", new_name))
		}
		names.out <- tax_table(x)[, pos.col]
		return(as.vector(names.out))
	}
}	

###Taxonomic cutoff function###
taxa_summary <- function(ps, taxrank, coff_val, KeepUnknown = TRUE,
	GroupByTaxon = TRUE, ...){
	valid.names <- validating.names(ps, taxrank = taxrank, ...)
	n.un <- which(colnames(tax_table(ps)) == taxrank)
	tax_table(ps)[, n.un] <- valid.names
	if(!KeepUnknown){
		tax_col <- as.vector(tax_table(ps)[, taxrank])
		keep_taxa <- tax_col != "Unknown"
		ps <- prune_taxa(keep_taxa, ps)
		valid.names <- as.vector(tax_table(ps)[, taxrank])
	}
	taxrel <- transform_sample_counts(ps, function(x) x / sum(x))
	df_asv <- as.data.frame(otu_table(taxrel))
	if(taxa_are_rows(taxrel)) df_asv <- t(df_asv)
	means.tmp <- colMeans(df_asv)
	df.tmp <- data.frame(
		ID = taxa_names(ps),
		Tax = valid.names,
		MeanRel = means.tmp,
		row.names = "ID",
		stringsAsFactors = F
	) %>%
	arrange(desc(MeanRel)) %>%
	mutate(CumSum = cumsum(MeanRel))
	if(GroupByTaxon){
		df.tmp <- df.tmp %>%
			group_by(Tax) %>%
			summarise(
				MeanRel = sum(MeanRel),
				ASV_Count = n()
			) %>%
			arrange(desc(MeanRel)) %>%
			mutate(CumSum = cumsum(MeanRel)) %>%
			ungroup() %>%
			mutate(RowLabel = paste0(ASV_Count, "_ASVs")) %>%
			as.data.frame()
		rownames(df.tmp) <- make.unique(df.tmp$RowLabel)
		df.tmp <- df.tmp[, c("Tax", "MeanRel", "CumSum")]
	}
	df.cut <- df.tmp %>% filter(CumSum <= coff_val)
	if(nrow(df.cut) < nrow(df.tmp)){
		df.cut <- bind_rows(df.cut, df.tmp[nrow(df.cut)+1, ])
	}
	top.tax <- df.cut$Tax
	tt <- tax_table(taxrel)
	idx <- as.character(tt[, taxrank])
	is.other <- !(idx %in% top.tax | idx == "Unknown")
	tt[is.other, taxrank] <- "Other"
	tax_table(taxrel) <- tt
	list(summary = df.cut, phy.lvl = taxrel)
}

###Optimal cutoff
optimal_cutoff <- function(start, by, ps, taxrank, ...){
	cutoffs <- seq(start, 1.0, by = by)
	screedata <- data.frame(
		Cutoff = cutoffs,
		NumTaxa = numeric(length(cutoffs)),
		CumAbundance = numeric(length(cutoffs))
	)
	for (i in seq_along(cutoffs)) {
  		tmp.ts <- taxa_summary(ps = ps, taxrank = taxrank,
			coff_val = cutoffs[i], ...
		)
  		screedata$NumTaxa[i] <- nrow(tmp.ts$summary)
		screedata$CumAbundance[i] <- max(tmp.ts$summary$CumSum)
	}
	tmp.scree <- screedata
	normal <- function(x){(x-min(x))/(max(x)-min(x))}
	tmp.scree$NumTaxa <- normal(screedata$NumTaxa)
	tmp.scree$CumAbundance <- normal(screedata$CumAbundance)
	tmp.scree$Dist <- abs(tmp.scree$NumTaxa - tmp.scree$CumAbundance)/sqrt(2)
	optCO <- tmp.scree$Cutoff[which.max(tmp.scree$Dist)]
	print(paste0(
		"Recommended cut-off after normalization and trade-off equality distance observation is: ",
		optCO
	))
	return(screedata)
}

#OC.phylum -> 0.95
#OC.class -> 0.95
#OC.order -> 0.80
#OC.family -> 0.75
#OC.genus -> 0.75

##Sample optimal_cutoff() and taxa_summary() outputs
OC.order <- optimal_cutoff(0.5, 0.05, pseq, "Order", KeepUnknown = F, NearestName = F, verbose = F)
OC.order[2:12,] <- OC.order[1:11,]; OC.order[1:11,] <- rep(0,3)
OC.order <- OC.order %>%
	mutate(slope = (Cutoff - lag(Cutoff))/(NumTaxa - lag(NumTaxa))) %>%
	mutate(slopeChange = slope/lag(slope))
OC.family <- optimal_cutoff(0.5, 0.05, pseq, "Family", KeepUnknown = F, NearestName = F, verbose = F)
OC.family[2:12,] <- OC.family[1:11,]; OC.family[1:11,] <- rep(0,3)
OC.family <- OC.family %>%
	mutate(slope = (Cutoff - lag(Cutoff))/(NumTaxa - lag(NumTaxa))) %>%
	mutate(slopeChange = slope/lag(slope))

NU.order <- taxa_summary(pseq, "Order", 0.80, KeepUnknown=F, NearestName = F, verbose = F)
NU.family <- taxa_summary(pseq, "Family", 0.80, KeepUnknown=F, NearestName = F, verbose = F)

##Plotting OC
ggplot(OC.order, aes(x = NumTaxa, y = CumAbundance))+
geom_line(color = "blue", lwd = 2)+
geom_point(shape = 17, size = 7)+
theme_minimal()

elbow.gg.order <- ggplot(OC.order, aes(x = NumTaxa, y = CumAbundance))+
geom_line(color = "blue", lwd = 2)+
theme_minimal()+
theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 20))+
geom_segment(
	aes(x = 24, xend = 24, y = 0, yend = 0.80),
	linetype = 'dashed',
	colour = 'red'
)+geom_segment(
	aes(x = 0, xend = 24, y = 0.80, yend = 0.80),
	linetype = 'dashed',
	colour = 'red'
)+geom_label(aes(label = NumTaxa))+
labs(x = "Number of Orders", y = "Cumulative Relative Abundance")

elbow.gg.family <- ggplot(OC.family, aes(x = NumTaxa, y = CumAbundance))+
geom_line(color = "blue", lwd = 2)+
theme_minimal()+
theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 20))+
geom_segment(
	aes(x = 29, xend = 29, y = 0, yend = 0.75),
	linetype = 'dashed',
	colour = 'red'
)+geom_segment(
	aes(x = 0, xend = 29, y = 0.75, yend = 0.75),
	linetype = 'dashed',
	colour = 'red'
)+geom_label(aes(label = NumTaxa))+
labs(x = "Number of Families", y = "Cumulative Relative Abundance")

jpeg("./SupplementaryFigure1.jpg", height = 4000, width = 4000, res = 300)
plot_grid(
	elbow.gg.order,
	elbow.gg.family,
	labels = c('(a)', '(b)'),
	ncol = 1
)
dev.off()

##Preparing ASV data
otu.genus <- as.data.frame(otu_table(NU.genus$phy.lvl))
colnames(otu.genus)[40] <- "sample-381"
otu.genus$Genus <- tax_table(NU.genus$phy.lvl)[,"Genus"]

##Preparing metadata
meta.genus <- as.data.frame(sample_data(NU.genus$phy.lvl))
rownames(meta.genus)[40] <- "sample-381"

##Collapsing data for hypothesis 1 (between sites)
h1_res.family <- list()
h1_sumlist.family <- list()
tmp.cp <- c(NU.family$summary$Tax)

for(i in tmp.cp){
	tmp.x <- otu.family[which(otu.family$Family==i),]
	tmp.Sum <- colSums(tmp.x[,-41])
	tmp.y <- data.frame(x = meta.family$species,
		y = tmp.Sum)
	h1_res.family[[i]] <- agricolae::kruskal(trt = tmp.y$x, y = tmp.y$y)
	hs <- h1_res.family[[i]]$statistics
	hs[6] <- i
	h1_sumlist.family[[i]] <- hs
}

size.df <- length(h1_res.family)
h1_family <- data.frame(
	name = character(size.df),
	chisq = numeric(size.o.df),
	p.val = numeric(size.df),
	sign = character(size.df)
)

for(i in 1:size.df){
	h1_family$name[i] = names(h1_res.family)[i]
	h1_family$chisq[i] = h1_res.family[[i]]$statistics$Chisq
	h1_family$p.val[i] = h1_res.family[[i]]$statistics$p.chisq
	if(h1_family$p.val[i] < 0.001){
		h1_family$sign[i] <- "***"
	} else if(h1_family$p.val[i] < 0.01){
		h1_family$sign[i] <- "**"
	} else if(h1_family$p.val[i] < 0.05){
		h1_family$sign[i] <- "*"
	} else if(h1_family$p.val[i] == 0.05){
		h1_family$sign[i] <- "."
	}
}

h1_family

##Collapsing data for hypothesis 2 (within sites)
h2_res.family <- list()
h2_sumlist.family <- list()
tmp.cp <- c(NU.family$summary$Tax)

for(i in tmp.cp){
	tmp.x <- otu.family[which(otu.family$Family==i),]
	tmp.Sum <- colSums(tmp.x[,-41])
	tmp.y <- data.frame(
		site = meta.family$species,
		type = meta.family$compartment,
		y = tmp.Sum
		)
	for(j in unique(tmp.y$site)){
		tmp.w <- subset(tmp.y, site == j)
		h2_res.family[[i]][[j]] <- agricolae::kruskal(
			trt = tmp.w$type,
			y = tmp.w$y
		)
		hs <- h2_res.family[[i]][[j]]$statistics
		hs[6] <- i
		h2_sumlist.family[[i]][[j]] <- hs
	}
}

size.df <- length(h2_res.family)
h2_family <- list()
for(k in names(h2_res.family[[1]])){
	h2_tmp <- data.frame(
		name = character(size.df),
		p.val = numeric(size.df),
		sign = character(size.df)
	)
	for(i in 1:size.df){
		h2_tmp$name[i] = names(h2_res.family)[i]
		h2_tmp$p.val[i] = h2_res.family[[i]][[k]]$statistics$p.chisq
		if(h2_tmp$p.val[i] < 0.001){
			h2_tmp$sign[i] <- "***"
		} else if(h2_tmp$p.val[i] < 0.01){
			h2_tmp$sign[i] <- "**"
		} else if(h2_tmp$p.val[i] < 0.05){
			h2_tmp$sign[i] <- "*"
		} else if(h2_tmp$p.val[i] == 0.05){
			h2_tmp$sign[i] <- "."
		}
	}
	h2_family[[k]] <- h2_tmp
}

h2_family

##Collapsing data for hypothesis 3 (among rhizospheres)
h3_res.family <- list()
h3_sumlist.family <- list()
tmp.cp <- c(NU.family$summary$Tax)

for(i in tmp.cp){
	tmp.x <- otu.family[which(otu.family$Family==i),]
	tmp.Sum <- colSums(tmp.x[,-41])
	tmp.y <- data.frame(
		site = meta.family$species,
		type = meta.family$compartment,
		y = tmp.Sum
		)
	tmp.w <- subset(tmp.y, type == "rhizosheath")
	h3_res.family[[i]] <- agricolae::kruskal(
		trt = tmp.w$site,
		y = tmp.w$y
	)
	hs <- h3_res.family[[i]]$statistics
	hs[6] <- i
	h3_sumlist.family[[i]] <- hs
}

size.df <- length(h3_res.family)
h3_family <- data.frame(
	name = character(size.df),
	chisq = numeric(size.df),
	p.val = numeric(size.df),
	sign = character(size.df)
)
for(i in 1:size.df){
	h3_family$name[i] = names(h3_res.family)[i]
	h3_fam	ily$chisq[i] = h3_res.family[[i]]$statistics$Chisq
	h3_family$p.val[i] = h3_res.family[[i]]$statistics$p.chisq
	if(h3_family$p.val[i] < 0.001){
		h3_family$sign[i] <- "***"
	} else if(h3_family$p.val[i] < 0.01){
		h3_family$sign[i] <- "**"
	} else if(h3_family$p.val[i] < 0.05){
		h3_family$sign[i] <- "*"
	} else if(h3_family$p.val[i] == 0.05){
		h3_family$sign[i] <- "."
	}
}

h3_family

##Collapsing data for hypothesis 4 (among bulk soils)
h4_res.family <- list()
h4_sumlist.family <- list()
tmp.cp <- c(NU.family$summary$Tax)

for(i in tmp.cp){
	tmp.x <- otu.family[which(otu.family$Family==i),]
	tmp.Sum <- colSums(tmp.x[,-41])
	tmp.y <- data.frame(
		site = meta.family$species,
		type = meta.family$compartment,
		y = tmp.Sum
		)
	tmp.w <- subset(tmp.y, type == "bulk_soil")
	h4_res.family[[i]] <- agricolae::kruskal(
		trt = tmp.w$site,
		y = tmp.w$y
	)
	hb <- h4_res.family[[i]]$statistics
	hb[6] <- i
	h4_sumlist.family[[i]] <- hb
}

size.df <- length(h4_res.family)
h4_family <- data.frame(
	name = character(size.df),
	chisq = numeric(size.df),
	p.val = numeric(size.df),
	sign = character(size.df)
)
for(i in 1:size.df){
	h4_family$name[i] = names(h4_res.family)[i]
	h4_family$chisq[i] = h4_res.family[[i]]$statistics$Chisq
	h4_family$p.val[i] = h4_res.family[[i]]$statistics$p.chisq
	if(h4_family$p.val[i] < 0.001){
		h4_family$sign[i] <- "***"
	} else if(h4_family$p.val[i] < 0.01){
		h4_family$sign[i] <- "**"
	} else if(h4_family$p.val[i] < 0.05){
		h4_family$sign[i] <- "*"
	} else if(h4_family$p.val[i] == 0.05){
		h4_family$sign[i] <- "."
	}
}

h4_family
c(mean(h4_family$chisq), mean(h4_family$p.val))


###Summary table of all four hypotheses median
tax_level <- c("Phylum", "Class", "Order", "Family", "Genus")
h1_all <- sapply(
	list(
		h1_phylum$p.val,
		h1_class$p.val,
		h1_order$p.val,
		h1_family$p.val,
		h1_genus$p.val
	), median
)
h2.1_Mcons <- sapply(
	list(
		h2_phylum$`M. conspicua`$p.val,
		h2_class$`M. conspicua`$p.val,
		h2_order$`M. conspicua`$p.val,
		h2_family$`M. conspicua`$p.val,
		h2_genus$`M. conspicua`$p.val
	), median
)
h2.2_Mhaag <- ... $same as `h2.1_Mcons`
h2.3_Mmeis <- ... $same as `h2.1_Mcons`
h2.4_Mlani <- ... $same as `h2.1_Mcons`
h3_MR <- ... #Same as `h1_all`
h4_BS <- ... #Same as `h1_all`
p_medians_bars <- data.frame(row.names = "tax_level", h1_all, h2.1_Mcons,
	h2.2_Mhaag, h2.3_Mmeis, h2.4_Mlani, h3_MR, h4_BS)

###Plotting this as a tileplot
tile_pmb <- ggplot(
	plot_pmb,
	aes(x = variable, y = Tax_Level, fill = value)
)+geom_tile(
	color = "white",
	linewidth = 0.8
)+geom_text(
	aes(label=formatC(
		value,
		format = 'e',
		digits = 1)),
	color = 'white', 
	size = 3
)+scale_fill_gradientn(
	colors = c("#FEE8C8", "#FC8D59", "#990000"),
	trans = "log10",
	breaks = 10^(-4:-1),
	labels = scales::scientific
)

tile_pmb

#levels(sample_data(NU.phylum$phy.lvl)$treatment)[4] <- "M. haageana MR" #Small correction

##Plotting barplots
sp.cols <- c("#F56C42", "#2995BC", "#43B385", "#FBB254")
library(ggh4x)
glom.phylum <- tax_glom(NU.phylum$phy.lvl, taxrank = "Phylum")
sample_names(NU.phylum$phy.lvl)[40] <- "sample-381"
bar.phylum <- plot_bar(
	glom.phylum,
	fill = "Phylum",
)+
geom_bar(aes(fill = Phylum),
	stat = 'identity',
	position = "stack"
)+
labs(y = "Relative frequency", x = "Soil Compartment")+
theme_minimal()+
theme(
	axis.title = element_text(size = 20),
	axis.text.x = element_text(size = 11, face = 'bold', angle = 90, hjust = 0),
	strip.text.x = element_text(size = 12, face = 'bold.italic'),
	legend.key.size = unit(8.6, 'mm'),
	legend.title = element_text(size = 16, face = 'bold')
)+facet_wrap2(
	~treatment,
	ncol = 8,
	scales = 'free',
	strip = strip_themed(
		background_x = elem_list_rect(
			fill = rep(sp.cols, each = 2)
		)
	), labeller = label_wrap_gen(width = 15)
)+theme(
	axis.text.y = element_blank(),
	axis.ticks.y = element_blank(),
	legend.key.width = unit(0.5, 'cm'),
	legend.text = element_text(size = 10)
)+scale_fill_manual(values = randomcoloR::distinctColorPalette(
	k = nrow(NU.phylum$summary)+1
	)
)

bar.phylum

###Plotting deltas
library("pheatmap")

mat_phylum <- as.data.frame(otu_table(glom.phylum))
mat_phylum$Tax <- as.vector(tax_table(glom.phylum)[,"Phylum"])
mat_phylum <- mat_phylum[mat_phylum$Tax != "Other",]
rownames(mat_phylum) <- mat_phylum$Tax
mat_phylum <- mat_phylum[,-41]
mat_phylum <- as.data.frame(t(mat_phylum))
mat_phylum <- mat_phylum %>% mutate(SampleID = rownames(.))
long_phylum <- mat_phylum %>% pivot_longer(
	cols = -SampleID,
	names_to = "phylum",
	values_to = "abundance"
)
meta.phylum$SampleID <- rownames(meta.phylum)
add_phylum <- as(meta.phylum, "data.frame") %>%
	select(
		SampleID,
		species,
		compartment) %>%
	rename(site = species)
long_phylum <- long_phylum %>%
	left_join(
		add_phylum,
		by = "SampleID"
	) %>% mutate(PairedID = substr(SampleID, 1, 9)) %>%
	select(
		PairedID,
		phylum,
		abundance,
		site,
		compartment) %>%
	pivot_wider(
		names_from = compartment,
		values_from = abundance,
		values_fill = 0
	) %>% mutate(delta = rhizosphere - bulk_soil)
delta_phylum <- long_phylum %>%
	select(
		phylum,
		PairedID,
		delta
	) %>%
	pivot_wider(
		names_from = PairedID, #Maybe change by site
		values_from = delta,
		values_fill = 0
	) %>%
	column_to_rownames("phylum") %>%
	as.matrix()
clust_phylum <- pheatmap(
	delta_phylum,
	cluster_rows = T,
	cluster_cols = T,
	main = "Δ frequency (rhizosphere − soil) per phylum"
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