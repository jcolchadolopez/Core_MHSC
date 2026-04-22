setwd("D:/CleanData")

library("phyloseq")
library("vegan")
library("ggplot2")
library("ggrepel")

#Loading the phyloseq object (abundance)
pseq <- readRDS("phyloseq_pseq.rds")

#Creating a matrix from otu_table
mat.otu <- t(otu_table(pseq)) #I need to t() as vegan uses samples as rows and OTUs as columns
class(mat.otu) <- "matrix" #This leads to a warning but is necessary to run rarefy()
df.rare <- rarecurve(mat.otu, steps = 1000, tidy = T)

#`tidy = T` creates a df that can be passed down into `ggplot2`-like functions

data.frame(
	"Site" = rownames(sample_data(pseq)),
	"trt" = sample_data(pseq)$treatment
)

ggplot(aes(
		x = Sample,
		y = Species,
		col = Site
	),
	data = df.rare)+
	geom_line()