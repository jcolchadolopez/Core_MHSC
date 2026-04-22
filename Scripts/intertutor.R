#### ESTE SCRIPT ES PARA EL TUTORAL ####
setwd("D:/CleanData")
library("microbiome")
library("microbiomeutilities")
library("phyloseq")
library("tidyverse")
library("vegan")
library("dplyr")

###Validating data function###
validating.names <- function(x, taxrank){
	pos.col <- which(colnames(tax_table(x))==taxrank)
	if(length(which(as.vector(tax_table(x)[,pos.col])=="uncultured")) == 0){
		print("All taxa have valid names")
		return(as.vector(tax_table(x)[,pos.col]))
	}
	else{
		error.pos <- which(as.vector(tax_table(x)[,pos.col])=="uncultured")
		print(
			paste0(
				"The following positions will be corrected: ",
				paste(error.pos, collapse = " , ")
			)
		)
		for(i in error.pos){
			print(paste0(i,": ", paste0(tax_table(x)[i,pos.col:2], collapse = "_")))
			tax_table(x)[i,pos.col] <- paste0(tax_table(x)[i,pos.col:2], collapse = "_")
		}
	new.names <- tax_table(x)[,pos.col]
	return(as.vector(new.names))
	}
}

###Class level - individual###
lvl.classis <- tax_glom(pseq, taxrank = "Classis")
taxa_names(lvl.classis) <- validating.names(lvl.classis, "Classis")
t10.classis <- names(sort(taxa_sums(lvl.classis), decreasing = T))[1:10]
taxa.classis <- tax_table(pseq)
taxa.classis[is.na(taxa.classis)] <- "Unknown"
taxa.classis <- taxa.classis[taxa.classis[,"Classis"] %in% t10.classis, ]
#sample_data(pseq)$sample.id <- rownames(sample_data(pseq))
#sample_data(pseq)$sample.id[40] <- "sample-381"
rseq.classis <- pseq
tax_table(rseq.classis) <- taxa.classis
colr.classis <- c("#ff8e00", "#8a5900", "#8e008e", "#008e00", "#ff0000",
	"#440044", "#00c0c0", "#6d4c00", "#4f3f00", "#650000")

rseq.classis <- transform_sample_counts(rseq.classis, function(x){ x/sum(x) }) 

bar.classis <- plot_bar(
	tax_glom(rseq.classis, taxrank = "Classis"),
	"compartment",
	fill = "Classis",
	facet_grid = ~species
)+
geom_bar(aes(fill = Classis),
	stat = 'identity',
	position = "stack"
)+
labs(y = "Relative frequency", x = "Soil Compartment")+
scale_fill_manual(values = colr.classis)+
theme(
	axis.title = element_text(size = 20),
	axis.text.y = element_text(size = 16, face = 'bold'),
	axis.text.x = element_text(size = 16, face = 'bold', angle = 0, hjust = 0.5),
	strip.text = element_text(size = 16, face = 'bold.italic'),
	legend.key.size = unit(8.6, 'mm'),
	legend.title = element_text(size = 16, face = 'bold'),
	legend.text = element_text(size = 12)
)+ggh4x::facet_wrap2(
	~species,
	ncol = 4,
	scales = 'free',
	strip = ggh4x::strip_themed(background_x = elem_list_rect(fill = sp.cols))
)

bar.classis

###Order level - individual ###
lvl.ordo <- tax_glom(rseq, taxrank = "Ordo")
taxa_names(lvl.ordo) <- validating.names(lvl.ordo, "Ordo")
t15.ordo <- names(sort(taxa_sums(lvl.ordo), decreasing = T))[1:15]
taxa.ordo <- tax_table(rseq)
taxa.ordo[is.na(taxa.ordo)] <- "Unknown"
taxa.ordo <- taxa.ordo[taxa.ordo[,"Ordo"] %in% t15.ordo, ]
rseq.ordo <- pseq
tax_table(rseq.ordo) <- taxa.ordo
colr.ordo <- c("#440044", "#8a5900", "#996000", "#413900", "#00c0c0",
	"#a76700", "#c47400", "#d37a00", "#e28100", "#990000", "#7c007c",
	"#6d4c00", "#4f3f00", "#767600", "#650000")

rseq.ordo <- transform_sample_counts(rseq.ordo, function(x){ x/sum(x) }) 

bar.ordo <- plot_bar(
	tax_glom(rseq.ordo, taxrank = "Ordo"),
	"compartment",
	fill = "Ordo",
	facet_grid = ~species
)+
geom_bar(aes(fill = Ordo),
	stat = 'identity',
	position = "stack"
)+
labs(y = "Relative frequency", x = "Soil Compartment")+
scale_fill_manual(values = colr.ordo)+
theme(
	axis.title = element_text(size = 20),
	axis.text.y = element_text(size = 16, face = 'bold'),
	axis.text.x = element_text(size = 16, face = 'bold', angle = 0, hjust = 0.5),
	strip.text = element_text(size = 16, face = 'bold.italic'),
	legend.key.size = unit(8.6, 'mm'),
	legend.title = element_text(size = 16, face = 'bold'),
	legend.text = element_text(size = 12)
)+ggh4x::facet_wrap2(
	~species,
	ncol = 4,
	scales = 'free',
	strip = ggh4x::strip_themed(background_x = elem_list_rect(fill = sp.cols))
)

bar.ordo

###Familia level - individual###
lvl.familia <- tax_glom(rseq, taxrank = "Familia")
tmp.names.familia <- validating.names(lvl.familia, "Familia")
tmp.names.familia[c(14,263)] <- c("Unknown_Gammaproteobacteria_Incertae_Sedis_Gammaproteobacteriia_Proteobacteria",
	"Unknown_Oxyphotobacteria_Incertae_Sedis_Cyanobacteriia_Cyanobacteria")
taxa_names(lvl.familia) <- tmp.names.familia
t15.familia <- names(sort(taxa_sums(lvl.familia), decreasing = T))[1:16]
taxa.familia <- tax_table(rseq)
taxa.familia[is.na(taxa.familia)] <- "Unknown"
taxa.familia <- taxa.familia[taxa.familia[,"Familia"] %in% t15.familia, ]
rseq.familia <- pseq
tax_table(rseq.familia) <- taxa.familia
colr.familia <- c("#4f3f00", "#008e00", "#7c007c", "#cc0000", "#00c0c0",
	"#c47400", "#d37a00", "#e28100", "#990000", "#6d4c00", "#5e4600",
	"#570057", "#767600", "#650000", "#690069")

rseq.familia <- transform_sample_counts(rseq.familia, function(x){ x/sum(x) }) 

bar.familia <- plot_bar(
	tax_glom(rseq.familia, taxrank = "Familia"),
	"compartment",
	fill = "Familia",
	facet_grid = ~species
)+
geom_bar(aes(fill = Familia),
	stat = 'identity',
	position = "stack"
)+
labs(y = "Relative frequency", x = "Soil Compartment")+
scale_fill_manual(values = colr.familia)+
theme(
	axis.title = element_text(size = 20),
	axis.text.y = element_text(size = 16, face = 'bold'),
	axis.text.x = element_text(size = 16, face = 'bold', angle = 0, hjust = 0.5),
	strip.text = element_text(size = 16, face = 'bold.italic'),
	legend.key.size = unit(8.6, 'mm'),
	legend.title = element_text(size = 16, face = 'bold'),
	legend.text = element_text(size = 12)
)+ggh4x::facet_wrap2(
	~species,
	ncol = 4,
	scales = 'free',
	strip = ggh4x::strip_themed(background_x = elem_list_rect(fill = sp.cols))
)

bar.familia

###Genus level - individual###
lvl.genus <- tax_glom(rseq, taxrank = "Genus")
taxa_names(lvl.genus) <- validating.names(lvl.genus, "Genus")
t20.genus <- names(sort(taxa_sums(lvl.genus), decreasing = T))[1:22]
taxa.genus <- tax_table(rseq)
taxa.genus[is.na(taxa.genus)] <- "Unknown"
taxa.genus <- taxa.genus[taxa.genus[,"Genus"] %in% t20.genus, ]
rseq.genus <- pseq
tax_table(rseq.genus) <- taxa.genus
colr.genus <- c("#ff0000", "#4f3f00", "#a76700", "#b66d00", "#5e4600",
	"#00c0c0", "#ffff00", "#bbbb00", "#7b5300", "#7c007c", "#8a5900",
	"#d37a00", "#e28100", "#990000", "#8e008e", "#400098", "#6d4c00",
	"#f08700", "#767600", "#650000")

rseq.genus <- transform_sample_counts(rseq.genus, function(x){ x/sum(x) }) 

bar.genus <- plot_bar(
	tax_glom(rseq.genus, taxrank = "Genus"),
	"compartment",
	fill = "Genus",
	facet_grid = ~species
)+
geom_bar(aes(fill = Genus),
	stat = 'identity',
	position = "stack"
)+
labs(y = "Relative frequency", x = "Soil Compartment")+
scale_fill_manual(values = colr.genus)+
guides(fill = guide_legend(ncol = 1))+
theme(
	axis.title = element_text(size = 20, face = 'bold'),
	axis.text.y = element_text(size = 16, face = 'bold'),
	axis.text.x = element_text(size = 16, face = 'bold', angle = 0, hjust = 0.5),
	strip.text = element_text(size = 20, face = 'bold.italic'),
	legend.key.size = unit(8.6, 'mm'),
	legend.title = element_text(size = 16, face = 'bold'),
	legend.text = element_text(size = 12)
)+
ggh4x::facet_wrap2(
	~species,
	ncol = 4,
	strip = ggh4x::strip_themed(background_x = elem_list_rect(fill = sp.cols)),
	scales = 'free'
)

bar.genus