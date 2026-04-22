setwd("D:/CleanData")
library("vegan")
library("dplyr")
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
#rownames(biom)[c(26,38)] <- c("sample-382", "sample-283")
biom <- biom[order(rownames(biom)),]

#DEPRECATED BECAUSE OF REMOVING 'NMR'
#meta[c(4,9,12,15,22,25,32,35,38,43,48,51),3:11] <- NA

#SUBSITUTE FOR:
#rm.nmr <- c(4,9,12,15,22,25,32,35,38,43,48,51)
#meta <- meta[-rm.nmr, ]
#biom <- biom[-rm.nmr, ]

meta$species <- as.factor(meta$species)
meta$compartment <- as.factor(meta$compartment)
taxa <- unique(biom)

set.seed(120)
####CHECAR SCRIPT: Adaptar para `pseq`!!!####
data.mds <- metaMDS(
	as(t(pseq@otu_table), "matrix"),
	k=5,
	distance = "bray",
	weakties = F,
	autotransform = F
)

stressplot(
	data.mds,
	main = paste(
		"Bray-Curtis (S=",
		round(
			data.mds$grstress,
			3
		),
		")",
		sep=""
	)
)

set.seed(120)
BC.envfit <- envfit(data.mds, sample_data(pseq), choices = c(1:5), permutations=9999, na.rm = T)
BC.spfit <- envfit(data.mds, as(t(pseq@otu_table), "matrix"), choices = c(1:5), permutations=999)

#Corrected, deprecate earliers
set.seed(120)
BC.envfit <- envfit(
	data.mds,
	pseq@sam_data,
	choices = c(1:5),
	permutations=9999,
	na.rm = T
)

set.seed(120)
BC.fmfit <- envfit(
	data.mds,
	as(
		t(tax_glom(pseq, "Family")@otu_table),
		"matrix"
	),
	choices = c(1:5),
	permutations=9999
)
rownames(BC.fmfit$vectors$arrows) <- tax_glom(pseq, "Family")@tax_table[,5]

set.seed(120)
BC.gnfit <- envfit(
	data.mds,
	as(
		t(tax_glom(pseq, "Genus")@otu_table),
		"matrix"
	),
	choices = c(1:5),
	permutations=999
)
rownames(BC.gnfit$vectors$arrows) <- tax_glom(pseq, "Genus")@tax_table[,6]

sample.scrs <- as.data.frame(scores(data.mds, display = "sites"))
sample.scrs <- cbind(
	sample.scrs,
	Host = pseq@sam_data$species,
	Zone = pseq@sam_data$compartment,
	Int = pseq@sam_data$treatment
)
sample.scrs

env.scores <- as.data.frame(
	scores(
		BC.envfit,
		display = "vectors"
	)
)
env.scores <- cbind(
	env.scores,
	env.variables = rownames(
		env.scores
	),
	pval = BC.envfit$vectors$pvals
)
sign.env.scrs <- subset(
	env.scores,
	pval<=0.001
)
sign.env.scrs

fam.scores <- as.data.frame(
	scores(
		BC.fmfit,
		display = "vectors"
	)
)
fam.scores <- cbind(
	fam.scores,
	fam.variables = rownames(
		fam.scores
	),
	pval = BC.fmfit$vectors$pvals
)
sign.fam.scrs <- subset(
	fam.scores,
	pval<=0.001
)
sign.fam.scrs

gen.scores <- as.data.frame(
	scores(
		BC.gnfit,
		display = "vectors"
	)
)
gen.scores <- cbind(
	gen.scores,
	gen.variables = rownames(
		gen.scores
	),
	pval = BC.gnfit$vectors$pvals
)
sign.gen.scrs <- subset(
	gen.scores,
	pval<=min(gen.scores$pval)
)
sign.gen.scrs

library("ggplot2")
library("ggrepel")
library("MetBrewer")

sp.cols <- met.brewer(name = "Egypt", n = 4, type = "discrete")

nmds.plot <- ggplot(
		sample.scrs,
		aes(
			x=NMDS1,
			y=NMDS4,
		)
	)+
	geom_point(
		aes(
			fill = as.factor(Host),
			shape = as.factor(Zone)
		),
		size = 4
	)+
	scale_fill_manual(
		name = "Host species (Locality)",
		values = sp.cols,
		labels = c(
			expression(bolditalic("M. conspicua")~bold("(Chazumba)")),
			expression(bolditalic("M. haageana")~bold("(Perote)")),
			expression(bolditalic("M. lanigera")~bold("(Nochixtlan)")),
			expression(bolditalic("M. meissneri")~bold("(Zapotitlan)"))
		)
	)+
	scale_shape_manual(
		name = "Soil Compartment",
		values = c(21, 24), #Rm = 22
		labels = c("Bulk Soil (BS)", "Rhizosheath (RS)")
	)+guides(
		fill = guide_legend(
			override.aes=list(shape=22),
			order = 1
		)#, 
#		shape = guide_legend(position = 'right')
	)+
	coord_fixed()+
	theme_bw()+
	theme(
		plot.title = element_text(
			face = "bold",
			size = 22,
			hjust = 0.5,
		), axis.title = element_text(
			face = "bold",
			size = 24
		), legend.text = element_text(
			size = 15,
			face = "italic"
		), legend.title = element_text(
			size = 16,
			face = "bold"
		), axis.text = element_text(
			size = 18
		), legend.position = "right"
	)+
	labs(
		title = "Bray-Curtis dissimilarity (S = 0.032)"
	)+
	geom_segment(
		data = sign.env.scrs,
		aes(
			x = 0,
			xend = NMDS1,
			y = 0,
			yend = NMDS4
		),
		arrow = arrow(
			length = unit(
				0.25,
				"cm")
			),
		colour = "black",
		lwd = 1.2
	)+
	ggrepel::geom_text_repel(
		data = sign.env.scrs,
		aes(
			x=NMDS1,
			y=NMDS4,
			label = c(
				"Humidity",
				"C mineralization",
				"pH",
				"Electrical conductivity",
				"Phosphomonoesterase",
				"B-Glucosidase",
				"Polyphenol oxidase",
#				"Chitinase",
				"Dehydrogenase"
			)
		), 
		colour = "red",
		fontface = "bold",
		cex = 5,
		direction = "both",
		segment.size = 0.25
	)
nmds.plot


#fill.legend <- ggplot(
#		sample.scrs,
#		aes(
#			x=NMDS1,
#			y=NMDS4,
#		)
#	)+
#	geom_point(
#		aes(
#			fill = as.factor(Host),
#			shape = as.factor(Zone)
#		),
#		size = 4
#	)+
#	scale_fill_manual(
#		name = "Host species (Locality)",
#		values = sp.cols,
#		labels = c("M. conspicua (Chazumba)", "M. haageana (Perote)",
#			"M. lanigera (Nochixtlan)", "M. meissneri (Zapotitlan)"
#		)
#	)+
#	scale_shape_manual(
#		name = "Soil Compartment",
#		values = c(21, 24), #Rm = 22
#		labels = c("Bulk Soil (BS)", "Rhizosheath (RS)")
#	)+guides(
#		fill = guide_legend(
#			override.aes=list(shape=22),
#			order = 1
#		),shape = 'none'
#	)+
#	coord_fixed()+
#	theme_bw()+
#	theme(
#		plot.title = element_text(
#			face = "bold",
#			size = 22,
#			hjust = 0.5,
#			),
#		axis.title = element_text(
#			face = "bold",
#			size = 24
#			),
#		legend.text = element_text(
#			size = 15,
#			face = "italic"
#		),
#		legend.title = element_text(
#			size = 16,
#			face = "bold"
#		), legend.justification = c(-1,-0.75),
#		axis.text = element_text(
#			size = 18
#		)
#	)
#fill.legend <- cowplot::get_legend(fill.legend)

nmds.env.plot <- nmds.plot+geom_segment(
	data = sign.env.scrs,
	aes(
		x = 0,
		xend = NMDS1,
		y = 0,
		yend = NMDS4
	),
	arrow = arrow(
		length = unit(
			0.25,
			"cm")
		),
	colour = "black",
	lwd = 1.2
)+ggrepel::geom_text_repel(
	data = sign.env.scrs,
	aes(
		x=NMDS1,
		y=NMDS4,
		label = env.variables
	), 
	colour = "red",
	fontface = "bold",
	cex = 5,
	direction = "both",
	segment.size = 0.25
)
nmds.env.plot

subset.fam.scrs <- subset(
	sign.fam.scrs,
	abs(NMDS4) >= 0.7
)
nmds.fam.plot <- nmds.plot+geom_segment(
	data = subset.fam.scrs,
	aes(
		x = 0,
		xend = NMDS1,
		y = 0,
		yend = NMDS4
	),
	arrow = arrow(
		length = unit(
			0.25,
			"cm")
		),
	colour = "black",
	lwd = 1.2
)+ggrepel::geom_text_repel(
	data = subset.fam.scrs,
	aes(
		x=NMDS1,
		y=NMDS4,
		label = fam.variables
	), 
	colour = "red",
	fontface = "bold",
	cex = 5,
	direction = "both",
	segment.size = 0.25
)
nmds.fam.plot

subset.gen.scrs <- subset(
	sign.gen.scrs,
	abs(NMDS4) >= 0.7
)
nmds.gen.plot <- nmds.plot+geom_segment(
	data = subset.gen.scrs,
	aes(
		x = 0,
		xend = NMDS1,
		y = 0,
		yend = NMDS4
	),
	arrow = arrow(
		length = unit(
			0.25,
			"cm")
		),
	colour = "blue",
	lwd = 1
)+ggrepel::geom_text_repel(
	data = subset.gen.scrs,
	aes(
		x=NMDS1,
		y=NMDS4,
		label = gen.variables
	), 
	colour = "purple",
	fontface = "bold",
	cex = 3,
	direction = "both",
	segment.size = 0.25
)
nmds.gen.plot

nmds.gen.plot+geom_segment(
	data = sign.env.scrs,
	aes(
		x = 0,
		xend = NMDS1,
		y = 0,
		yend = NMDS4
	),
	arrow = arrow(
		length = unit(
			0.25,
			"cm")
		),
	colour = "black",
	lwd = 1
)+ggrepel::geom_text_repel(
	data = sign.env.scrs,
	aes(
		x=NMDS1,
		y=NMDS4,
		label = env.variables
	), 
	colour = "red",
	fontface = "bold",
	cex = 3,
	direction = "both",
	segment.size = 0.25
)

#####PCA
Newest.PCA <- ggplot(
		df.pca,
		aes(
			x=Dim.1,
			y=Dim.3,
		)
	)+
	geom_point(
		aes(
			fill = as.factor(Species),
			shape = as.factor(Zone)
		),
		size = 4
	)+
	scale_fill_manual(
		values = c(
			"blue",
			"green",
			"purple",
			"orange"
		)
	)+
	scale_shape_manual(
		values = c(21, 22)
	)+
	coord_fixed()+
	theme_bw()+
	theme(
		plot.title = element_text(
			face = "bold",
			size = 22,
			hjust = 0.5,
			),
		axis.title = element_text(
			face = "bold",
			size = 24
			),
		legend.position = "bottom",
		legend.text = element_text(
			size = 15,
			face = "italic"
		),
		legend.title = element_text(
			size = 16,
			face = "bold"
		),
		axis.text = element_text(
			size = 18
		)
	)+
	labs(
		fill = "Host species",
		shape = "Soil region",
		x = "PC1 (51.3%)",
		y = "PC3 (10.0%)"		
	)+
	guides(
		fill = guide_legend(
			override.aes=list(
				shape=21
			)
		)
	)+
	geom_segment(
		data = flechas,
		aes(
			x = 0,
			xend = Dim.1,
			y = 0,
			yend = Dim.3
		),
		arrow = arrow(
			length = unit(
				0.25,
				"cm")
			),
		colour = "black",
		lwd = 1.2
	)+
	ggrepel::geom_text_repel(
		data = flechas,
		aes(
			x=Dim.1,
			y=Dim.3,
			label = c(
				"Hum.",
				"CO2",
				"pH",
				"EC",
				"PME",
				"B-Glu",
				"POx",
				"Chi.",
				"DH"
			)
		), 
		colour = "red",
		fontface = "bold",
		cex = 4,
		direction = "both",
		segment.size = 0.25
	)
New.PCA
Newest.PCA

plot_grid(New.PCA, Newest.PCA, nrow=2)


##BOXPLOTS##
bpdb <- data.frame(Sample = rownames(sample.scrs), Host = sample.scrs$Host, Zone = sample.scrs$Zone, sample.scrs[,1:5])
bpdb[,1] <- sample.scrs$Host
set.seed(120)
N1.bp <- ggplot(bpdb, aes(x=Zone, y=NMDS1, fill=Host))+
	geom_boxplot()+
	geom_jitter()+
	scale_fill_manual(
		values = c("green", "blue", "orange", "purple"))+
	facet_wrap(.~Host, nrow=1)+
	theme(axis.title.x = element_blank(),
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		legend.position = "none"
		)

N2.bp <- ggplot(bpdb, aes(x=Zone, y=NMDS2, fill=Host))+
	geom_boxplot()+
	geom_jitter()+
	scale_fill_manual(
		values = c("green", "blue", "orange", "purple"))+
	facet_wrap(.~Host, nrow=1)+
	theme(axis.title.x = element_blank(),
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		legend.position = "none"
		)

N3.bp <- ggplot(bpdb, aes(x=Zone, y=NMDS3, fill=Host))+
	geom_boxplot()+
	geom_jitter()+
	scale_fill_manual(
		values = c("green", "blue", "orange", "purple"))+
	facet_wrap(.~Host, nrow=1)+
	theme(axis.title.x = element_blank(),
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		legend.position = "none"
		)

N4.bp <- ggplot(bpdb, aes(x=Zone, y=NMDS4, fill=Host))+
	geom_boxplot()+
	geom_jitter()+
	scale_fill_manual(
		values = c("green", "blue", "orange", "purple"))+
	facet_wrap(.~Host, nrow=1)+
	theme(axis.title.x = element_blank(),
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		legend.position = "none"
		)

N5.bp <- ggplot(bpdb, aes(x=Zone, y=NMDS5, fill=Host))+
	geom_boxplot()+
	geom_jitter()+
	scale_fill_manual(
		values = c("green", "blue", "orange", "purple"))+
	facet_wrap(.~Host, nrow=1)+
	scale_x_discrete(labels = c("BS", "SC", "RZ"))+
	theme(legend.position = "none")


#library("cowplot")
plot_grid(N1.bp, N2.bp, N3.bp, N4.bp, N5.bp,
	nrow = 5
)

 

library("cowplot")
plot_grid(nmds.MC, nmds.MH, nmds.MM, nmds.ML,
	labels = c("M. conspicua", "M. haageana", "M. meissneri", "M. lanigera"),
	label_size = 12
)

library("reshape2")
phyla <- read.csv(
	"L2_rels.csv",
	header=T
)

phyla$X.OTU.ID <- c("Haageana BS", "Haageana RZ",
	"Haageana BS", "Haageana NS", "Haageana RZ",
	"Haageana BS", "Haageana RZ",
	"Haageana BS", "Haageana NS", "Haageana RZ",
	"Haageana BS", "Haageana NS",
	"Conspicua NS", "Conspicua NS", "Conspicua NS",
	"Meissneri BS", "Meissneri RZ",
	"Meissneri BS", "Meissneri RZ",
	"Meissneri BS", "Meissneri NS", "Meissneri RZ",
	"Meissneri BS", "Meissneri NS", "Meissneri RZ",
	"Meissneri BS", "Meissneri RZ",
	"Lanigera BS", "Lanigera RZ",
	"Lanigera BS", "Lanigera NS", "Lanigera RZ",
	"Lanigera BS", "Lanigera RZ",
	"Lanigera BS", "Lanigera NS", "Lanigera RZ",
	"Lanigera BS", "Lanigera NS", "Lanigera RZ",
	"")

library("plyr")

phean <- ddply(
	phyla,
	"X.OTU.ID",
	numcolwise(mean),
	na.rm = T
)

phylt <- melt(
	phean[-1,],
	na.rm = T
)

colnames(phylt) <- c("Sample", "Phylum", "value")

phylt$Sample <- factor(
	phylt$Sample,
	levels = unique(phylt$Sample)
)

barphylt <- ggplot(
		phylt,
		aes(
			x = Sample,
			fill = Phylum,
			y = value
		)
	)+geom_bar(
		stat = "identity",
		colour = "black"
	)+theme(
		axis.text.x = element_text(
			angle = 90,
			size = 14,
			colour = "black",
			vjust = 0.5,
			hjust = 1,
			face= "bold"
		), 
		axis.title.y = element_text(
			size = 16,
			face = "bold"
		),
		legend.title = element_text(
			size = 16,
			face = "bold"
		),
		legend.text = element_text(
			size = 12,
			face = "bold",
			colour = "black"
		),
		axis.text.y = element_text(
			colour = "black",
			size = 12,
			face = "bold"
		)
	)+scale_y_continuous(expand = c(0,0))+
	labs(x = "",
		y = "Relative Abundance",
		fill = "Phylum"
	)

barphylt

###Kruskal-Wallis###
KW.phyla <- melt(
	phyla[-1,],
	na.rm = T
)

colnames(KW.phyla) <- c("Sample", "Phylum", "value")

List.phyla <- levels(KW.phyla$Phylum)

KW.lists <- c()
KW.listing <- function(n){
	KW.ind <- KW.phyla[KW.phyla$Phylum==List.phyla[n],]
	KW.res <- kruskal(
		y = KW.ind$value,
		trt = KW.ind$Sample,
		alpha = 0.05,
		p.adj = "bonferroni",
		group = T
	)
	KW.st <- KW.res$statistics
	KW.st$sign. <- ifelse(KW.st$p.chisq < 0.001,
		"***",
		ifelse(KW.st$p.chisq < 0.01,
			"**",
			ifelse(KW.st$p.chisq < 0.05,
				"*",
				"n.s."
			)
		)
	)
	colnames(KW.st)[4] <- paste(
		"s.",
		List.phyla[n],
		sep = ''
	)
	KW.lt <- KW.res$groups
	colnames(KW.lt)[1] <- List.phyla[n]
	KW.sub <- c()
	KW.sub[[1]] <- KW.st
	KW.sub[[2]] <- KW.lt
	return(KW.sub)
}
	
set.seed(120)
as_mat <- as.matrix(biom)
Compa <- anosim(
	as_mat,
	as.factor(meta$compartment),
	distance = "bray",
	permutations = 999
)
Compa

Morra <- anosim(
	as_mat,
	as.factor(meta$species),
	distance = "bray",
	permutations = 999
)
Morra

meta$sp_x_cp <- paste(meta$species,
	"*",
	meta$compartment,
	sep=''
)
Sola <- anosim(
	as_mat,
	as.factor(meta$sp_x_cp),
	distance = "bray",
	permutations = 999
)


library("eulerr")
#library("BiocManager")
#BiocManager::install("microbiome")
#devtools::install_github("microsud/microbiomeutilities")
library("microbiome")
library("devtools")
library("microbiomeutilities")
library("dplyr")

####Haageana####
Haa.meta <- meta[meta$species=="haageana",]
veemo <- read.table(
	

Haa.seq <- phyloseq(
	otu_table(
		biom,
		taxa_are_rows = F
	), sample_data(Haa.meta)
)

Haa.rel <- microbiome::transform(
	Haa.seq,
	"compositional"
)

Haa.list_venn <- unique(as.character(meta(Haa.rel)$compartment))
Haa.list_venn

Haa.list_core <- c()
for (n in Haa.list_venn){
	ps.sub <- subset_samples(
		Haa.rel,
		compartment == n
	)
	core_m <- core_members(
		ps.sub,
		detection = 0.01,
		prevalence = 0.75
	)
	print(
		paste0(
			"No. of core taxa in ",
			n,
			" : ",
			length(core_m)
		)
	)
	Haa.list_core[[n]] <- core_m
}

mycols.C <- c(
	bulk_soil = "green",
	rhizosphere = "red",
	near_soil = "blue"
)

####Conspicua####
Coi.meta <- meta[meta$species=="conspicua",]

Coi.seq <- phyloseq(
	otu_table(
		biom,
		taxa_are_rows = F
	), sample_data(Coi.meta)
)

Coi.rel <- microbiome::transform(
	Coi.seq,
	"compositional"
)

Coi.list_venn <- unique(as.character(meta(Coi.rel)$compartment))
Coi.list_venn

Coi.list_core <- c()
for (n in Coi.list_venn){
	ps.sub <- subset_samples(
		Coi.rel,
		compartment == n
	)
	core_m <- core_members(
		ps.sub,
		detection = 0.01,
		prevalence = 0.75
	)
	print(
		paste0(
			"No. of core taxa in ",
			n,
			" : ",
			length(core_m)
		)
	)
	Coi.list_core[[n]] <- core_m
}

####Lanigera####
Lai.meta <- meta[meta$species=="lanigera",]

Lai.seq <- phyloseq(
	otu_table(
		biom,
		taxa_are_rows = F
	), sample_data(Lai.meta)
)

Lai.rel <- microbiome::transform(
	Lai.seq,
	"compositional"
)

Lai.list_venn <- unique(as.character(meta(Lai.rel)$compartment))
Lai.list_venn

Lai.list_core <- c()
for (n in Lai.list_venn){
	ps.sub <- subset_samples(
		Lai.rel,
		compartment == n
	)
	core_m <- core_members(
		ps.sub,
		detection = 0.01,
		prevalence = 0.75
	)
	print(
		paste0(
			"No. of core taxa in ",
			n,
			" : ",
			length(core_m)
		)
	)
	Lai.list_core[[n]] <- core_m
}

HVP <- plot(venn(Haa.list_core), fills = mycols.C)
CVP <- plot(venn(Coi.list_core), fills = mycols.C)
LVP <- plot(venn(Lai.list_core), fills = mycols.C)

####Compartimentos####
Rhz.meta <- meta[meta$species!="meissneri",]
Rhz.meta <- meta[meta$compartment=="rhizosphere",]

Rhz.seq <- phyloseq(
	otu_table(
		biom,
		taxa_are_rows = F
	), sample_data(Rhz.meta)
)

Rhz.rel <- microbiome::transform(
	Rhz.seq,
	"compositional"
)

Rhz.list_venn <- unique(as.character(meta(Rhz.rel)$species))
Rhz.list_venn

Rhz.list_core <- c()
for (n in Rhz.list_venn){
	ps.sub <- subset_samples(
		Rhz.rel,
		species == n
	)
	core_m <- core_members(
		ps.sub,
		detection = 0.001,
		prevalence = 0.75
	)
	print(
		paste0(
			"No. of core taxa in ",
			n,
			" : ",
			length(core_m)
		)
	)
	Rhz.list_core[[n]] <- core_m
}

mycols.S <- c(
	haageana = "blue",
	conspicua = "green",
	lanigera = "purple"
)

RVP <- plot(venn(Rhz.list_core), fills = mycols.S)

####Haageana####
Haa.list_total <- c()
for (n in Haa.list_venn){
	ts.sub <- subset_samples(
		Haa.rel,
		compartment == n
	)
	core_t <- core_members(
		ts.sub,
		detection = 1/10000000,
		prevalence = 1/100,
		include.lowest = T
	)	
	print(
		paste0(
			"No. of core taxa in ",
			n,
			" : ",
			length(core_t)
		)
	)
	Haa.list_total[[n]] <- core_t
}

tHVP <- plot(venn(Haa.list_total), fills = mycols.C)

####Conspicua####
Coi.list_total <- c()
for (n in Coi.list_venn){
	ts.sub <- subset_samples(
		Coi.rel,
		compartment == n
	)
	core_t <- core_members(
		ts.sub,
		detection = 1/10000000,
		prevalence = 1/100,
		include.lowest = T
	)	
	print(
		paste0(
			"No. of core taxa in ",
			n,
			" : ",
			length(core_t)
		)
	)
	Coi.list_total[[n]] <- core_t
}

tCVP <- plot(venn(Coi.list_total), fills = mycols.C)

####Lanigera####
Lai.list_total <- c()
for (n in Lai.list_venn){
	ts.sub <- subset_samples(
		Lai.rel,
		compartment == n
	)
	core_t <- core_members(
		ts.sub,
		detection = 1/10000000,
		prevalence = 1/100,
		include.lowest = T
	)	
	print(
		paste0(
			"No. of core taxa in ",
			n,
			" : ",
			length(core_t)
		)
	)
	Lai.list_total[[n]] <- core_t
}

tLVP <- plot(venn(Lai.list_total), fills = mycols.C)

####Compartimentos####
Rhz.list_total <- c()
for (n in Rhz.list_venn){
	ts.sub <- subset_samples(
		Rhz.rel,
		species == n
	)
	core_t <- core_members(
		ts.sub,
		detection = 1/10000000,
		prevalence = 1/100,
		include.lowest = T
	)	
	print(
		paste0(
			"No. of core taxa in ",
			n,
			" : ",
			length(core_t)
		)
	)
	Rhz.list_total[[n]] <- core_t
}

tRVP <- plot(venn(Rhz.list_total), fills = mycols.S)
