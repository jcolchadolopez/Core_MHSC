setwd("D:/CleanData")
library("microbiome")
library("microbiomeutilities")
library("phyloseq")
library("tidyverse")
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

rownames(meta)[29] <- "sample-341"
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

#rm.nmr <- c(4,9,12,15,22,25,32,35,38,43,48,51)
#meta <- meta[-rm.nmr, ]
#biom <- biom[-rm.nmr, ]

meta$species <- as.factor(meta$species)
meta$compartment <- as.factor(meta$compartment)
meta$treatment <- as.factor(paste0(meta$species, "_", meta$compartment))
levels(meta$treatment) <- c("M. conspicua BS", "M. conspicua MR",
					"M. haageana BS", "M. haagean MR",
					"M. lanigera BS", "M. lanigera MR",
					"M. meissneri BS", "M. meissneri MR"
	)
taxa <- unique(biom)

set.seed(120)

samples_df <- meta
head(samples_df)

otu_mat <- as.matrix(t(biom))
head(otu_mat)

taxa <- read.csv(
	"TAXAQ2.csv",
	header=T,
	row.names = 1
)

tbiom <- as.data.frame(t(biom))
tbiom$ROW <- rownames(tbiom)
taxa$ROW <- rownames(taxa)
taxa_filt <- merge(taxa, tbiom, by = "ROW")
rownames(taxa_filt) <- taxa_filt$ROW
taxa_filt <- taxa_filt[order(rownames(taxa_filt)),]
taxa_filt <- taxa_filt[,-c(3:ncol(taxa_filt))]
tax_mat <- taxa_filt %>%
	tidyr::separate(Taxon,
		into = c(
			"Dominium", "Phylum", "Classis",
			"Ordo", "Familia", "Genus", "Species"
		),
		sep = ';'
	) %>%	dplyr::mutate(across("Dominium",
		~ifelse(
			. == "",
			"NA",
			substr(., 4, nchar(.))
		)
	)) %>% dplyr::mutate(
		across(c(
			"Phylum", "Classis", "Ordo",
			"Familia", "Genus", "Species"
		),
		~ifelse(
			. == "",
			"NA",
			substr(., 5, nchar(.))
		)
	))
tax_mat <- as.matrix(tax_mat[,-1])

OTU <- otu_table(otu_mat, taxa_are_rows = T)
TAX <- tax_table(tax_mat)
samples <- sample_data(samples_df)

pseq <- phyloseq(OTU, TAX, samples)

write

sample_names(pseq)
rank_names(pseq)
sample_variables(pseq)

#Saving this object
base::saveRDS(pseq, "./phyloseq_pseq.rds")

bseq <- prune_taxa(taxa_sums(pseq) > 1, pseq)
rseq <- merge_samples(bseq, "treatment")
rseq <- transform_sample_counts(rseq, function(x){ x/sum(x) })
sample_data(rseq)[,2] <- rep(c("BS", "MR"), 4)
sample_data(rseq)[,1] <- c(
	rep("M. conspicua",2),
	rep("M. haageana",2),
	rep("M. lanigera",2),
	rep("M. meissneri",2)
)

#Saving this one too
base::saveRDS(rseq, "./phyloseq_rseq.rds")

###INICIO DEL SALTO###
bar.D <- plot_bar(
	rseq,
	"compartment",
	fill = "Dominium",
	facet_grid = ~host
)+
geom_bar(aes(color=Dominium, fill = Dominium),
	stat = 'identity',
	position = "stack")+
labs(y = "Relative frequency", x = "Soil Compartment")
###FIN DEL SALTO###

phylos <- tax_glom(rseq, taxrank = "Phylum")
taxa_names(phylos) <- tax_table(phylos)[,"Phylum"]
t10.phylos <- names(sort(taxa_sums(phylos), decreasing = T))[1:10]
taxa.phylos <- tax_table(rseq)
taxa.phylos[is.na(taxa.phylos)] <- "Unknown"
taxa.phylos[!taxa.phylos [,"Phylum"] %in% t10.phylos & taxa.phylos[,"Phylum"] != "Unknown", "Phylum"] <- "Other"
rseq.phylos <- rseq
tax_table(rseq.phylos) <- taxa.phylos
bar.P <- plot_bar(
	tax_glom(rseq.phylos, taxrank = "Phylum"),
	"compartment",
	fill = "Phylum",
	facet_grid = ~species
)+
geom_bar(aes(color=Phylum, fill = Phylum),
	stat = 'identity',
	position = "stack"
)+
labs(y = "Relative frequency", x = "Soil Compartment")+
scale_fill_manual(values = c(rainbow(n = 11), "gray"))

classe <- tax_glom(rseq, taxrank = "Classis")
taxa_names(classe) <- c(tax_table(classe)[,"Classis"][1:37], "uncultured_Desulfobacterota",
				tax_table(classe)[,"Classis"][39:45], "uncultured_Armatimonadota",
				tax_table(classe)[,"Classis"][47:length(taxa_names(classe))])
t10.classe <- names(sort(taxa_sums(classe), decreasing = T))[1:10]
taxa.classe <- tax_table(rseq)
taxa.classe[is.na(taxa.classe)] <- "Unknown"
taxa.classe[!taxa.classe[,"Classis"] %in% t10.classe & taxa.classe[,"Classis"] != "Unknown", "Classis"] <- "Other"
rseq.classe <- rseq
tax_table(rseq.classe) <- taxa.classe
bar.C <- plot_bar(
	tax_glom(rseq.classe, taxrank = "Classis"),
	"compartment",
	fill = "Classis",
	facet_grid = ~species
)+
geom_bar(aes(fill = Classis),
	stat = 'identity',
	position = "stack"
)+
labs(y = "Relative frequency", x = "Soil Compartment")+
scale_fill_manual(values = c(rainbow(n = 10)[1:7], "gray", rainbow(n = 10)[8:10], "black"))+
theme(
	axis.title = element_text(size = 20),
	axis.text = element_text(size = 16),
	strip.text = element_text(size = 16, face = 'italic')
)

ordem <- tax_glom(rseq, taxrank = "Ordo")
temp.taxa <- function(x, rank){
	x.rnk <- which(colnames(tax_table(x)) == rank)
	tmp.D <- sum(duplicated(tax_table(x)[,x.rnk]))
	if(tmp.D > 0){
		tmp.which <- which(tax_table(x)[,x.rnk] == "uncultured")
		for(i in tmp.which){
			tax_table(x)[,x.rnk][i] <- ifelse(tax_table(x)[i,x.rnk-1] == "uncultured",
				paste0("uncultured_prov.", i),
				paste0(
					tax_table(x)[i,x.rnk],
					"_",
					tax_table(x)[i,x.rnk-1]
				)
			)
		}
	}
	return(tax_table(x)[,x.rnk])
}
taxa_names(ordem) <- temp.taxa(ordem, "Ordo")
t15.ordem <- names(sort(taxa_sums(ordem), decreasing = T))[1:15]
taxa.ordem <- tax_table(rseq)
taxa.ordem[is.na(taxa.ordem)] <- "Unknown"
taxa.ordem[!taxa.ordem[,"Ordo"] %in% t15.ordem & taxa.ordem[,"Ordo"] != "Unknown", "Ordo"] <- "Other"
rseq.ordem <- rseq
tax_table(rseq.ordem) <- taxa.ordem
colres <- sample(rainbow(n = 15))

bar.O <- plot_bar(
	tax_glom(rseq.ordem, taxrank = "Ordo"),
	"compartment",
	fill = "Ordo",
	facet_grid = ~host
)+
geom_bar(aes(fill = Ordo),
	stat = 'identity',
	position = "stack"
)+
labs(y = "Relative frequency", x = "Soil Compartment")+
scale_fill_manual(values = c(colres[1:7], "gray", colres[8:14], "black", colres[15]))+
theme(
	axis.title = element_text(size = 20, face = 'bold'),
	axis.text.y = element_text(size = 16, face = 'bold', ),
	axis.text.x = element_text(size = 16, face = 'bold', angle = 0, hjust = 0.5),
	strip.text = element_text(size = 20, face = 'bold.italic'),
	legend.key.size = unit(8.6, 'mm'),
	legend.title = element_text(size = 16, face = 'bold'),
	legend.text = element_text(size = 12)
)+
ggh4x::facet_wrap2(
	~host,
	ncol = 4,
	strip = strip_themed(background_x = elem_list_rect(fill = c("green", "blue", "purple", "orange"))),
	
)
bar.O

####Rarefactions
library("vegan")
read.count <- data.table(as(sample_data(bseq), 'data.frame'),
	TotalReads = sample_sums(bseq),
	keep.rownames = T)
setnames(read.count, 'rn', 'SampleID')
read.depth <- ggplot(read.count, aes(TotalReads))+
	geom_histogram()+
	ggtitle("Sequencing Depth")+
	theme(title = element_text(size = 24),
		axis.title = element_text(size = 20),
		axis.text = element_text(size = 16)
	)
vegan::rarecurve(as.data.frame(t(otu_table(bseq))), step = 1000, label = T)


###VENN
library("eulerr")
trts <- rownames(meta(rseq))
list_omne <- c()
list_over <- c()
for(n in 1:nrow(meta(rseq))){
	ps.sub <- subset_samples(rseq, rownames(otu_table(rseq)) %in% trts[n])
	ps.sub <- prune_taxa(colSums(otu_table(ps.sub)) > 0, ps.sub)
	omne.m <- colnames(otu_table(ps.sub))
	list_over[[1]][n] <- print(paste0("No. of core taxa in ", trts[n], " : ", length(omne.m)))
	list_omne[[n]] <- omne.m
	list_comp.1 <- list(list_over, list_omne)
}
plot(venn(list_omne))

#Compartimentos
comp.trts <- unique(meta(rseq)$compartment)
clst_omne <- c()
for(n in comp.trts){
	c.sub <- subset_samples(rseq, compartment == n)
	c.sub <- prune_taxa(colSums(otu_table(c.sub)) > 0, c.sub)
	omne.c <- colnames(otu_table(c.sub))
	print(paste0("No. of core taxa in ", n, " : ", length(omne.c)))
	clst_omne[[n]] <- omne.c
}
p.omne <- plot(venn(clst_omne))

########################
#Combinado: F > 0
host.trts <- unique(meta(rseq)$species)
comp.trts <- unique(meta(rseq)$compartment)
list_omne <- c()
for(n in host.trts){
		m.sub <- subset_samples(rseq, species == n)
		m.sub <- prune_taxa(colSums(otu_table(m.sub)) > 0, m.sub)
		m.sub <- subset_samples(m.sub, compartment == "MR")
		m.sub <- prune_taxa(colSums(otu_table(m.sub)) > 0, m.sub)
		omne.m <- colnames(otu_table(m.sub))
		name <- paste0(n, "_", "MR")
		print(paste0("No. of core taxa in ", name, " : ", length(omne.m)))
		list_omne[[name]] <- omne.m
}
sp.cols <- MetBrewer::met.brewer(name = "Egypt", n = 4, type = "discrete")
names(list_omne) <- c("M. conspicua", "M. haageana", "M. lanigera", "M. meissneri")
p.list <- ggvenn(list_omne, fill_color = sp.cols, show_percentage = F)+
	labs(title = "MR shared core taxa among Mammillaria species")+
	theme(plot.title = element_text(
			size = 30,
			face = "bold",
			hjust = 0.5
			),
		legend.text = element_text(face = "italic")
		)
p.list
core_asvs <- Reduce(intersect, list_omne)
taxa_core <- tax_table(rseq)[core_asvs,]
6 Acidobacteria; 14 Actinobacteriota, 1 Bacteroidota, 7 Chloroflexi, 1 Cyanobacteria,
1 Entotheonellaeota, 1 Firmicutes, 2 Gemmatimonadota, 2 Myxococca, 1 Nitrospirota,
1 Planctomycetota, 14 Proteobacteria

45 Proteobacteria:Alphaproteobacteria:Rhizobiales
22 


install.packages(
	pkgs = "C:/Users/RosasLab Acer/Downloads/Tax4Fun2_1.1.5.tar.gz",
	repos = NULL,
	source = TRUE
)
library("Tax4Fun2")
str(taxa_core)
t4f_obj <- Tax4Fun(as.data.frame(tax_table), reference " 





#Combinado
host.trts <- unique(meta(rseq)$species)
comp.trts <- unique(meta(rseq)$compartment)
list_omne <- c()
for(n in host.trts){
		m.sub <- subset_samples(rseq, species == n)
		m.sub <- prune_taxa(colSums(otu_table(m.sub)) > 0.001, m.sub)
#	for(m in comp.trts){
		m.sub <- subset_samples(m.sub, compartment == "MR")
		m.sub <- prune_taxa(colSums(otu_table(m.sub)) > 0.001, m.sub)
		omne.m <- colnames(otu_table(m.sub))
		name <- paste0(n, "_MR")
		print(paste0("No. of core taxa in ", name, " : ", length(omne.m)))
		list_omne[[name]] <- omne.m
#	}
}

sp.cols <- MetBrewer::met.brewer(name = "Egypt", n = 4, type = "discrete")
names(list_omne) <- c("M. conspicua", "M. haageana", "M. lanigera", "M. meissneri")
q.list <- ggvenn(list_omne, fill_color = sp.cols)+
	labs(title = "MR shared core taxa among Mammillaria species")+
	theme(plot.title = element_text(
			size = 30,
			face = "bold",
			hjust = 0.5
			),
		legend.text = element_text(face = "italic")
		)
q.list

core_asvs <- Reduce(intersect, list_omne)
taxa_core <- tax_table(rseq)[core_asvs,]

c("e21846d75cda56b0a6d7f2ef8b52e890", "5ed84bd03b448df039cd8a6bc73012a7", "dc76f594eef06fecb1810c695e1bfb3c")
################################33



#Individual: cons
comp.trts <- unique(meta(rseq)$compartment)
list_lani <- c()

l.sub <- subset_samples(rseq, species == "lanigera")
l.sub <- prune_taxa(colSums(otu_table(l.sub)) > 0.01, l.sub)
l.sub <- subset_samples(l.sub, compartment == comp.trts[3])
l.sub <- prune_taxa(colSums(otu_table(l.sub)) > 0.01, l.sub)
omne.l <- colnames(otu_table(l.sub))
list_lani[[comp.trts[3]]] <- omne.l

plot(venn(list_cons), main = "conspicua", fill = c("#66BB11", "#80FF66", "#80FF00"))
dev.new()
plot(venn(list_haag), main = "haageana", fill = c("#663C91", "#8080E6", "#808080"))
dev.new()
plot(venn(list_meis), main = "meissneri", fill = c("#E68E11", "#FFD266", "#FFD200"))
plot(venn(list_meis), main = "meissneri", fill = c("#A63C51", "#C080A6", "#C08040"))
dev.new()
cons.vd <- venn.diagram(
	list_cons,
	filename = "cons_vd.png",
	print.mode = c('raw', 'percent'),
	category.names = c("BS", "NMR", "MR"),
	fill = c("#66BB11", "#80FF66", "#80FF00")
)

haag.vd <- venn.diagram(
	list_haag,
	filename = "haag_vd.png",
	print.mode = c('raw', 'percent'),
	category.names = c("BS", "NMR", "MR"),
	fill = c("#663C91", "#8080E6", "#808080")
)

lani.vd <- venn.diagram(
	list_lani,
	filename = "lani_vd.png",
	print.mode = c('raw', 'percent'),
	category.names = c("BS", "NMR", "MR"),
	fill = c("#E68E11", "#FFD266", "#FFD200")
)

meis.vd <- venn.diagram(
	list_meis,
	filename = "meis_vd.png",
	print.mode = c('raw', 'percent'),
	category.names = c("BS", "NMR", "MR"),
	fill = c("#A63C51", "#C080A6", "#C08040")
)

library("UpSetR")
colores <- c("#000000", "#FF0000", "#FF00FF", "#800080", "#800000", "#FF0080", "#808080")
upset(
	as.data.frame(fromList(list_spp)),
	nsets = 12,
	main.bar.color = colores,
	sets.bar.color = rep(c("green", "blue", "purple", "orange"), 3),
	queries = list(list(
		query = intersects,
		params = list('conspicua_BS', 'conspicua_NS', 'conspicua_RH'),
		color = "green",
		active = T
	))
)