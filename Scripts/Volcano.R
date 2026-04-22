library(phyloseq)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(tidyverse)

# Transform into DSEq object
dseq <- phyloseq_to_deseq2(pseq, ~ species + compartment)

# Estimate size factors
dseq <- estimateSizeFactors(dseq, type = "poscounts")

# Differential analysis
dseq <- DESeq(dseq, test = "Wald", fitType = "parametric")
dseq.res <- results(
	dseq,
	contrast = c("compartment", "rhizosheath", "bulk_soil"),
	alpha = 0.05)

# Prepare the data for the volcano plot
res.df <- as.data.frame(dseq.res) %>%
	rownames_to_column("ASV") %>%
	left_join(as.data.frame(
		as.data.frame(tax_table(pseq)) %>% rownames_to_column("ASV")
		)
	) %>%
	mutate(
		log2FoldChange = ifelse(
			is.na(log2FoldChange),
			0,
			log2FoldChange
		),
		padj = ifelse(
			is.na(padj),
			1,
			padj
		),
		Significance = case_when(
			padj < 0.05 & log2FoldChange > 1 ~ "Up in RH",
			padj < 0.05 & log2FoldChange < -1 ~ "Up in BS",
			TRUE ~ "Non-significant"
		),
		Label = ifelse(
			padj < 0.001 & abs(log2FoldChange) > 2,
			coalesce(Genus, ASV),
			NA
		)
	)

write.csv(res.df, "./Volcano_data.csv")


# 4. Volcano plot
Don_Goyo <- ggplot(
	res.df,
	aes(
		x = log2FoldChange, y = -log10(padj),
		color = Significance,
		label = Label
	))+geom_point(alpha = 0.6, size = 2.5)+
	geom_vline(
		xintercept = c(-1, 1),
		linetype = "dashed",
		alpha = 0.5
	)+geom_hline(
		yintercept = -log10(0.05),
		linetype = "dashed",
		alpha = 0.5
	)+scale_color_manual(
		values = c("Up in RH" = "red",
			"Up in BS" = "blue",
			"Non-significant" = "gray"
		)
	)+geom_text_repel(
		max.overlaps = 20,
		size = 3,
		box.padding = 0.5,
		segment.color = "grey50"
	)+labs(
		x = "log2(Fold Change) Rhizosheath (RH) vs Bulk Soil (BS)",
		y = "-log10(Adjusted p-value)",
		title = "Differential Abundance: Compartment RH vs BS",
#		subtitle = "Controlled for species effects (n=4 species)"
	)+theme_minimal(base_size = 12)+
	theme(
		legend.position = "top",
		panel.grid.minor = element_blank(),
		plot.title = element_text(face = "bold")
	)+scale_shape_manual(
			values = c("Up in RH" = 24,
			"Up in BS" = 25,
			"Non-significant" = 21
		)
	)

Don_Goyo