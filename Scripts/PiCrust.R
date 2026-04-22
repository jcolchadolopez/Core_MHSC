#Putative functional annotation based on taxonomical identities

#Packages to load:
library("tidyverse")
library("phyloseq")
library("rstatix")
library("vegan")
library("picante")
library("kableExtra")
library("reticulate")
#BiocManager::install("ALDEx2")
library("ALDEx2")
library("data.table")

setwd("D:/CleanData")
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

#Creating phyloseq object
#1) Creating the OTU Table
otu_mat <- as.matrix(t(biom))
# head(otu_mat)
OTU <- otu_table(t(otu_mat), taxa_are_rows = T)

#2) Creating the Sample Data
samples <- sample_data(meta)

#3) Creating the Taxonomy Table:
taxa <- read.csv(
	"TAXAQ2.csv",
	header=T,
	row.names = 1
)

tbiom <- as.data.frame(biom)
tbiom$ROW <- rownames(tbiom)
taxa$ROW <- rownames(taxa)
taxa_filt <- merge(taxa, tbiom, by = "ROW")
rownames(taxa_filt) <- taxa_filt$ROW
taxa_filt <- taxa_filt[order(rownames(taxa_filt)),]
taxa_filt <- taxa_filt[,-c(3:ncol(taxa_filt))]
tax_mat <- taxa_filt %>%
	separate(Taxon,
		into = c(
			"Dominium", "Phylum", "Classis",
			"Ordo", "Familia", "Genus", "Species"
		),
		sep = ';'
	) %>%	mutate(across("Dominium",
		~ifelse(
			. == "",
			"NA",
			substr(., 4, nchar(.))
		)
	)) %>% mutate(
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
TAX <- tax_table(tax_mat)

#4) Creating the Phylogenetic Tree:
tree <- read.tree("tree.nwk")

#5) Creating the DNAStringSet:
refseqs <- Biostrings::readDNAStringSet("refseqs.fasta")
DNA <- refseq(refseqs)

#6) Joining them all into a phyloseq class object
pseq <- phyloseq(OTU, TAX, samples, DNA, tree)

#Cleaning the `phyloseq` class object
# WARNING: There were still 3 ARCHAEA sequences!
# Must check again my cleaning process in Qiime2 ###### URGENT #####
for(i in colnames(TAX)[-1]){
	ps.clean <- subset_taxa(pseq, Dominium == "Bacteria") %>%
		subset_taxa(!i %in% c("Chloroplast")) %>%
		subset_taxa(!i %in% c("mitochondria"))
}
ps.clean

# Filtering out singletons
# LOG: Further 247 ASVs were removed
ps.filt <- filter_taxa(ps.clean, function(x){sum(x > 0) >= 1 }, prune = T)
ps.filt

#Preparing the data for PICRUSt
dir.create("picrust")
#---- phyloseq_counts_to_df
phyloseq_counts_to_df <- function(pso){
	df.0 <- as.data.frame(otu_table(pso))
	df.0 <- cbind("#OTU ID" = rownames(df.0), df.0)
	return(df.0)
}
#---- write_counts_to_file
write_counts_to_file <- function(pso, file){
	df.1 <- phyloseq_counts_to_df(pso)
	write.table(df,
              file=file,
              sep = "\t",
              row.names = F
	)
}
#---- write_seqs_to_file
write_seqs_to_file <- function(pso, file){
	df.2 <- as.data.frame(refseq(pso))
	unlink(file)
	name2 <- rownames(df.2)
	for(i in nrow(df.2)){
		cat(
			paste0(">", name2[i]),
			file=file,
			sep="\n",
			append=TRUE
		)
		cat(
			df.2[i,1],
			file=file,
			sep="\n",
			append=TRUE
		)
	}
	message(
		paste0(
			"Wrote ",
			nrow(df.2),
			" ASV sequences to file ",
			file
		)
	)
}
write_counts_to_file(ps.filt, file = "picrust/raw_counts.tsv")
write_seqs_to_file(ps.filt, file = "picrust/seqs.fna")

#Installing PICRUSt2
system('mkdir -p picrust2 ; cd picrust2 ; rm -Rf results; picrust2_pipeline.py -s seqs.fna -i raw_counts.tsv -p 4 -o results --stratified --verbose')