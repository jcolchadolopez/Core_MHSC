require("vegan")
require("ggplot2")
require("tidyverse")

paths.table <- read_tsv("./picrust_out/pathways_out/path_abun_unstrat.tsv")

paths.table <- paths.table[, sort(names(paths.table))]

sample.info <- data.frame(
	Sample = colnames(paths.table)[-1],
	species = as.vector(pseq@sam_data$species),
	compartment = as.vector(pseq@sam_data$compartment),
	index = substr(colnames(paths.table)[-1], 8, 9)
)

#To add metadata to the PiCrust2 output
paths.table <- paths.table %>%
	pivot_longer(
		cols = -pathway,
		names_to = "Sample",
		values_to = "Abundance"
	) %>%
	left_join(
		sample.info,
		by = "Sample"
	) %>%  #This is to control the number of decimals printed
	mutate(
		across(
			where(is.numeric),
			~ num(., digits = 6)
		)
	)

#Now, we group by average and by sum per pathway to find the elbow points
paths.sum <- paths.table %>%
	group_by(pathway) %>%
	summarise(
		mean_abund = mean(Abundance, na.rm = T),
		sum_abund = sum(Abundance, na.rm = T)
	)

#Then, first sort by average to add rank
paths.sum.means <- paths.sum %>%
	arrange(desc(mean_abund)) %>%
	mutate(
		rank_mean = row_number(),
		cum_mean_prop = cumsum(mean_abund/sum(mean_abund))
	)

#Now we rank by sum
paths.sum.agg <- paths.sum %>%
	arrange(desc(sum_abund)) %>%
	mutate(
		rank_sum = row_number(),
		cum_sum_prop = cumsum(sum_abund/sum(sum_abund))
	)

#Add ranks to the main tibble
paths.sum <- paths.sum %>%
	left_join(
		select(paths.sum.means, "pathway", "rank_mean"),
		by = "pathway"
	) %>%
	left_join(
		select(paths.sum.agg, "pathway", "rank_sum"),
		by = "pathway"
	) %>%
	mutate(
		ave_rank = (rank_mean + rank_sum)/2
	) %>%
	arrange(desc(ave_rank)) %>%
	select(!c(rank_mean, rank_sum))

#Now, let's calculate the elbow point
#install.packages("smerc")
elbow.points.est <- as.data.frame(
		smerc::elbow_point(paths.sum$ave_rank, paths.sum$mean_abund),
		ncol = 3,
		row.names = "Mean"
	)
elbow.points.est <- rbind(
	elbow.points.est,
	as.data.frame(
		smerc::elbow_point(paths.sum$ave_rank, paths.sum$sum_abund),
		ncol = 3,
		row.names = "Sum"
	)
)
elbow.points.est[,1] <- c(
	paths.sum$pathway[elbow.points.est[1,2]],
	paths.sum$pathway[elbow.points.est[2,2]]	
)
elbow.points.est

#And plot it
elbow.paths <- ggplot()+
	geom_line(
		data = paths.sum.means,
		aes(x = rank_mean, y = cum_mean_prop),
		linewidth = 2, col = 'blue'
	)+
	geom_line(
		data = paths.sum.agg,
		aes(x = rank_sum, y = cum_sum_prop),
		linewidth = 1, col = 'red'
	)+
	labs(
		title = "Pathways Abundance Elbow Point",
		x = "Pathway Rank (descending order)",
		y = "Cummulative Abundance (red - sum, blue - mean)"
	)+
	geom_hline(
		yintercept = paths.sum.means$cum_mean_prop[37],
		col = "black",
		linetype = "dotted"
	)+
	annotate(
		geom = "text",
		x = 50,
		y = 0.25,
		label = "Elbow point at rank 37 (PWY-5686)",
		col = "black",
		hjust = 0,
		vjust = 0
	)
elbow.paths
		 
#Adding rank to the pathways table
paths.table <- paths.table %>%
	left_join(
		select(paths.sum, pathway, ave_rank),
		by = "pathway"
	) %>%
	arrange(-desc(ave_rank))


##################
##Let's redo it for each species by separate
#First 'M. haageana'
paths.haag.sum <- paths.table %>%
	filter(species == "M. haageana") %>%
	group_by(pathway) %>%
	summarise(
		mean_abund = mean(Abundance, na.rm = T),
		sum_abund = sum(Abundance, na.rm = T)
	) %>%	arrange(desc(mean_abund)) %>%
	mutate(
		rank_mean = row_number(),
		cum_mean_prop = cumsum(mean_abund/sum(mean_abund))
	) %>% arrange(desc(sum_abund)) %>%
	mutate(
		rank_sum = row_number(),
		cum_sum_prop = cumsum(mean_abund/sum(mean_abund))
	) %>%
	mutate(
		ave_rank = (rank_mean + rank_sum)/2
	) %>%
	arrange(-desc(ave_rank))
	
#Second 'M. conspicua'
paths.cons.sum <- paths.table %>%
	filter(species == "M. conspicua") %>%
	group_by(pathway) %>%
	summarise(
		mean_abund = mean(Abundance, na.rm = T),
		sum_abund = sum(Abundance, na.rm = T)
	) %>%	arrange(desc(mean_abund)) %>%
	mutate(
		rank_mean = row_number(),
		cum_mean_prop = cumsum(mean_abund/sum(mean_abund))
	) %>% arrange(desc(sum_abund)) %>%
	mutate(
		rank_sum = row_number(),
		cum_sum_prop = cumsum(mean_abund/sum(mean_abund))
	) %>%
	mutate(
		ave_rank = (rank_mean + rank_sum)/2
	) %>%
	arrange(-desc(ave_rank))

#Third comes 'M. meissneri'
paths.meis.sum <- paths.table %>%
	filter(species == "M. meissneri") %>%
	group_by(pathway) %>%
	summarise(
		mean_abund = mean(Abundance, na.rm = T),
		sum_abund = sum(Abundance, na.rm = T)
	) %>%	arrange(desc(mean_abund)) %>%
	mutate(
		rank_mean = row_number(),
		cum_mean_prop = cumsum(mean_abund/sum(mean_abund))
	) %>% arrange(desc(sum_abund)) %>%
	mutate(
		rank_sum = row_number(),
		cum_sum_prop = cumsum(mean_abund/sum(mean_abund))
	) %>%
	mutate(
		ave_rank = (rank_mean + rank_sum)/2
	) %>%
	arrange(-desc(ave_rank))

#Lastly, 'M. lanigera'
paths.lani.sum <- paths.table %>%
	filter(species == "M. lanigera") %>%
	group_by(pathway) %>%
	summarise(
		mean_abund = mean(Abundance, na.rm = T),
		sum_abund = sum(Abundance, na.rm = T)
	) %>%	arrange(desc(mean_abund)) %>%
	mutate(
		rank_mean = row_number(),
		cum_mean_prop = cumsum(mean_abund/sum(mean_abund))
	) %>% arrange(desc(sum_abund)) %>%
	mutate(
		rank_sum = row_number(),
		cum_sum_prop = cumsum(mean_abund/sum(mean_abund))
	) %>%
	mutate(
		ave_rank = (rank_mean + rank_sum)/2
	) %>%
	arrange(-desc(ave_rank))

## Looking for the smerc::elbow_point, these where the results
# M. haageana: at rank 35 (both methods)
# M. conspicua: at rank 35 (both methods)
# M. meissneri: at rank 35 (both methods)
# M. lanigera: at rank 37 (both methods)

#Extract the name of metabolic pathways per species
unique(paths.haag.sum$pathway)[1:35] -> haag.labs
unique(paths.cons.sum$pathway)[1:35] -> cons.labs
unique(paths.meis.sum$pathway)[1:35] -> meis.labs
unique(paths.lani.sum$pathway)[1:37] -> lani.labs #Up to 37 according to elbow_point()
unique(paths.table)[1:37] -> path.labs #This is for global

#Join them in a single list
list.labs <- list(
	M.haageana = haag.labs,
	M.conspicua = cons.labs,
	M.meissneri = meis.labs,
	M.lanigera = lani.labs,
	MHSC = path.labs
)

#Create a binary matrix
matrix.labs.names <- sort(unique(c(
	haag.labs,
	cons.labs,
	meis.labs,
	lani.labs,
	path.labs
)))
matrix.vals.haag <- which(matrix.labs.names %in% haag.labs)
matrix.vals.cons <- which(matrix.labs.names %in% cons.labs)
matrix.vals.meis <- which(matrix.labs.names %in% meis.labs)
matrix.vals.lani <- which(matrix.labs.names %in% lani.labs)
matrix.vals.MHSC <- which(matrix.labs.names %in% path.labs)
binary.venn <- matrix(nrow = 40, ncol = 5)
rownames(binary.venn) <- matrix.labs.names
colnames(binary.venn) <- names(list.labs)
matrix.vals.haag
binary.venn[matrix.vals.haag,1] <- 1; binary.venn[-matrix.vals.haag,1] <- 0
binary.venn[matrix.vals.cons,2] <- 1; binary.venn[-matrix.vals.cons,2] <- 0
binary.venn[matrix.vals.meis,3] <- 1; binary.venn[-matrix.vals.meis,3] <- 0
binary.venn[matrix.vals.lani,4] <- 1; binary.venn[-matrix.vals.lani,4] <- 0
binary.venn[matrix.vals.MHSC,5] <- 1; binary.venn[-matrix.vals.MHSC,5] <- 0
diff.venn <- c()
for(i in 1:nrow(binary.venn)){
	if(sum(binary.venn[i,])<5){diff.venn[i] <- i}
	else{diff.venn[i] <- NA}
}
binary.venn[!is.na(diff.venn),]

##################
##Let's redo it for each species by separate compartments
#First 'M. haageana'
paths.haag.bs <- paths.table %>%
	filter(species == "M. haageana", compartment == "bulk_soil") %>%
	group_by(pathway) %>%
	summarise(
		mean_abund = mean(Abundance, na.rm = T),
		sum_abund = sum(Abundance, na.rm = T)
	) %>%	arrange(desc(mean_abund)) %>%
	mutate(
		rank_mean = row_number(),
		cum_mean_prop = cumsum(mean_abund/sum(mean_abund))
	) %>% arrange(desc(sum_abund)) %>%
	mutate(
		rank_sum = row_number(),
		cum_sum_prop = cumsum(mean_abund/sum(mean_abund))
	) %>%
	mutate(
		ave_rank = (rank_mean + rank_sum)/2
	) %>%
	arrange(-desc(ave_rank))

paths.haag.rs <- paths.table %>%
	filter(species == "M. haageana", compartment == "rhizosheath") %>%
	group_by(pathway) %>%
	summarise(
		mean_abund = mean(Abundance, na.rm = T),
		sum_abund = sum(Abundance, na.rm = T)
	) %>%	arrange(desc(mean_abund)) %>%
	mutate(
		rank_mean = row_number(),
		cum_mean_prop = cumsum(mean_abund/sum(mean_abund))
	) %>% arrange(desc(sum_abund)) %>%
	mutate(
		rank_sum = row_number(),
		cum_sum_prop = cumsum(mean_abund/sum(mean_abund))
	) %>%
	mutate(
		ave_rank = (rank_mean + rank_sum)/2
	) %>%
	arrange(-desc(ave_rank))
	
#Second 'M. conspicua'
paths.cons.bs <- paths.table %>%
	filter(species == "M. conspicua", compartment == "bulk_soil") %>%
	group_by(pathway) %>%
	summarise(
		mean_abund = mean(Abundance, na.rm = T),
		sum_abund = sum(Abundance, na.rm = T)
	) %>%	arrange(desc(mean_abund)) %>%
	mutate(
		rank_mean = row_number(),
		cum_mean_prop = cumsum(mean_abund/sum(mean_abund))
	) %>% arrange(desc(sum_abund)) %>%
	mutate(
		rank_sum = row_number(),
		cum_sum_prop = cumsum(mean_abund/sum(mean_abund))
	) %>%
	mutate(
		ave_rank = (rank_mean + rank_sum)/2
	) %>%
	arrange(-desc(ave_rank))

paths.cons.rs <- paths.table %>%
	filter(species == "M. conspicua", compartment == "rhizosheath") %>%
	group_by(pathway) %>%
	summarise(
		mean_abund = mean(Abundance, na.rm = T),
		sum_abund = sum(Abundance, na.rm = T)
	) %>%	arrange(desc(mean_abund)) %>%
	mutate(
		rank_mean = row_number(),
		cum_mean_prop = cumsum(mean_abund/sum(mean_abund))
	) %>% arrange(desc(sum_abund)) %>%
	mutate(
		rank_sum = row_number(),
		cum_sum_prop = cumsum(mean_abund/sum(mean_abund))
	) %>%
	mutate(
		ave_rank = (rank_mean + rank_sum)/2
	) %>%
	arrange(-desc(ave_rank))

#Third comes 'M. meissneri'
paths.meis.bs <- paths.table %>%
	filter(species == "M. meissneri", compartment == "bulk_soil") %>%
	group_by(pathway) %>%
	summarise(
		mean_abund = mean(Abundance, na.rm = T),
		sum_abund = sum(Abundance, na.rm = T)
	) %>%	arrange(desc(mean_abund)) %>%
	mutate(
		rank_mean = row_number(),
		cum_mean_prop = cumsum(mean_abund/sum(mean_abund))
	) %>% arrange(desc(sum_abund)) %>%
	mutate(
		rank_sum = row_number(),
		cum_sum_prop = cumsum(mean_abund/sum(mean_abund))
	) %>%
	mutate(
		ave_rank = (rank_mean + rank_sum)/2
	) %>%
	arrange(-desc(ave_rank))

paths.meis.rs <- paths.table %>%
	filter(species == "M. meissneri", compartment == "rhizosheath") %>%
	group_by(pathway) %>%
	summarise(
		mean_abund = mean(Abundance, na.rm = T),
		sum_abund = sum(Abundance, na.rm = T)
	) %>%	arrange(desc(mean_abund)) %>%
	mutate(
		rank_mean = row_number(),
		cum_mean_prop = cumsum(mean_abund/sum(mean_abund))
	) %>% arrange(desc(sum_abund)) %>%
	mutate(
		rank_sum = row_number(),
		cum_sum_prop = cumsum(mean_abund/sum(mean_abund))
	) %>%
	mutate(
		ave_rank = (rank_mean + rank_sum)/2
	) %>%
	arrange(-desc(ave_rank))

#Then the fourth, 'M. lanigera'
paths.lani.bs <- paths.table %>%
	filter(species == "M. lanigera", compartment == "bulk_soil") %>%
	group_by(pathway) %>%
	summarise(
		mean_abund = mean(Abundance, na.rm = T),
		sum_abund = sum(Abundance, na.rm = T)
	) %>%	arrange(desc(mean_abund)) %>%
	mutate(
		rank_mean = row_number(),
		cum_mean_prop = cumsum(mean_abund/sum(mean_abund))
	) %>% arrange(desc(sum_abund)) %>%
	mutate(
		rank_sum = row_number(),
		cum_sum_prop = cumsum(mean_abund/sum(mean_abund))
	) %>%
	mutate(
		ave_rank = (rank_mean + rank_sum)/2
	) %>%
	arrange(-desc(ave_rank))

paths.lani.rs <- paths.table %>%
	filter(species == "M. lanigera", compartment == "rhizosheath") %>%
	group_by(pathway) %>%
	summarise(
		mean_abund = mean(Abundance, na.rm = T),
		sum_abund = sum(Abundance, na.rm = T)
	) %>%	arrange(desc(mean_abund)) %>%
	mutate(
		rank_mean = row_number(),
		cum_mean_prop = cumsum(mean_abund/sum(mean_abund))
	) %>% arrange(desc(sum_abund)) %>%
	mutate(
		rank_sum = row_number(),
		cum_sum_prop = cumsum(mean_abund/sum(mean_abund))
	) %>%
	mutate(
		ave_rank = (rank_mean + rank_sum)/2
	) %>%
	arrange(-desc(ave_rank))

#Lastly, the overall 'Mammillaria haageana Species Complex' or 'MHSC'
paths.MHSC.bs <- paths.table %>%
	filter(compartment == "bulk_soil") %>%
	group_by(pathway) %>%
	summarise(
		mean_abund = mean(Abundance, na.rm = T),
		sum_abund = sum(Abundance, na.rm = T)
	) %>%	arrange(desc(mean_abund)) %>%
	mutate(
		rank_mean = row_number(),
		cum_mean_prop = cumsum(mean_abund/sum(mean_abund))
	) %>% arrange(desc(sum_abund)) %>%
	mutate(
		rank_sum = row_number(),
		cum_sum_prop = cumsum(mean_abund/sum(mean_abund))
	) %>%
	mutate(
		ave_rank = (rank_mean + rank_sum)/2
	) %>%
	arrange(-desc(ave_rank))

paths.MHSC.rs <- paths.table %>%
	filter(compartment == "rhizosheath") %>%
	group_by(pathway) %>%
	summarise(
		mean_abund = mean(Abundance, na.rm = T),
		sum_abund = sum(Abundance, na.rm = T)
	) %>%	arrange(desc(mean_abund)) %>%
	mutate(
		rank_mean = row_number(),
		cum_mean_prop = cumsum(mean_abund/sum(mean_abund))
	) %>% arrange(desc(sum_abund)) %>%
	mutate(
		rank_sum = row_number(),
		cum_sum_prop = cumsum(mean_abund/sum(mean_abund))
	) %>%
	mutate(
		ave_rank = (rank_mean + rank_sum)/2
	) %>%
	arrange(-desc(ave_rank))

## Looking for the smerc::elbow_point, these where the results
# M. haageana - BS: at rank 204 (both methods) - N25: 36 - N50: 85
# M. haageana - RS: at rank 219(both methods) - N25: 39 - N50: 93
# Differences of +3, +8
# M. conspicua - BS: at rank 203 (both methods) - N25: 36 - N50: 84
# M. conspicua - RS: at rank 210 (both methods) - N25: 38 - N50: 89
# Differences of +2, +5
# M. meissneri - BS: at rank 207  (both methods) - N25: 37 - N50: 87
# M. meissneri - RS: at rank 210 (both methods) - N25: 38 - N50: 91
# Differences of +1, +4
# M. lanigera - BS: at rank 215 (both methods) - N25: 38 - N50: 91
# M. lanigera - RS: at rank 210 (both methods) - N25: 39 - N50: 94
# Differences of +1, +3
# MHSC - BS: at rank 211 (both methods) - N25: 37 - N50: 87 - 2.312184e8
# MHSC - RS: at rank 214 (both methods) - N25: 39 - N50: 92 - 2.475864e8
# Differences of +2, +5

# Are the 40 overall found important at local level too?
# Compare the ranks they occupy
# In general, both lists coincide on their first 34 positions
# BUT there are notable exceptions
which(unique(paths.haag.bs$pathway) %in% matrix.labs.names)
#Enriched rare metabolism: 21, 28, 34 & 35
# 21: PWY-7323 - O-Antigen; biotic stress (virus, protists) 
# 28: PWY-7357 - Thiamine salvage; C syntrophy
# 34: PWY-7790 - UMP; damaged nucleic acids 
# 35: PWY-7791 - UMP; damaged nucleic acids
which(unique(paths.haag.rs$pathway) %in% matrix.labs.names)
which(unique(paths.cons.bs$pathway) %in% matrix.labs.names)
# 34: PWY-7790 - UMP; damaged nucleic acids 
# 35: PWY-7791 - UMP; damaged nucleic acids
which(unique(paths.cons.rs$pathway) %in% matrix.labs.names)
which(unique(paths.meis.bs$pathway) %in% matrix.labs.names)
which(unique(paths.meis.rs$pathway) %in% matrix.labs.names)
# 34: PWY-7345: sucrose degradation; copiotrophic metabolism
which(unique(paths.lani.bs$pathway) %in% matrix.labs.names)
which(unique(paths.lani.rs$pathway) %in% matrix.labs.names)
which(unique(paths.MHSC.bs$pathway) %in% matrix.labs.names)
which(unique(paths.MHSC.bs$pathway) %in% matrix.labs.names)


######################BARPLOTS#########################
top.paths <- paths.table %>%
	filter(pathway %in% matrix.labs.names) %>%
	group_by(pathway, species, compartment) %>%
	summarise(
		mean_abund = mean(Abundance, na.rm = T),
		sd_abund = sd(Abundance, na.rm = T),
		.groups = 'drop'
	) %>%
	arrange(mean_abund) %>%
	mutate(
		across(
			where(is.numeric),
			~ num(., digits = 6)
		)
	) %>%
	mutate(pathway = factor(pathway, levels = unique(pathway)))

library(ggh4x)
jpeg("./Heatmaps/barplot_top_pathways.jpeg", height = 4000, width = 4000, res = 300)
top.plot <- ggplot(top.paths,
		aes(x = mean_abund, y = pathway, fill = compartment)
	)+geom_col(position = position_dodge(width = 0.8))+
	geom_errorbar(
		aes(
			xmin = mean_abund - sd_abund,
			xmax = mean_abund + sd_abund
		),
		width = 0.3,
		position = position_dodge(width = 0.8)
	)+facet_wrap2(
		~species,
		scales = 'free_y',
		strip = strip_themed(
			background_x = elem_list_rect(fill = sp.cols)
		), labeller = label_wrap_gen(width = 13)
	)+
	labs(
		x = "Predicted abundance",
		y = NULL,
		fill = "Compartment"
	)+
	theme_bw(base_size = 14)+
	theme(strip.text.x = element_text(size = 12, face = 'bold.italic'))+
	scale_fill_manual(
		values = c(
			"bulk_soil" = 'sienna',
			"rhizosheath" = 'wheat'
		), labels = c(
			"bulk_soil" = "Bulk Soil",
			"rhizosheath" = "Rhizosheath"
		)
	)
top.plot
dev.off()

######################DELTA#########################
delta.paths <- paths.table %>%
	pivot_wider(
		names_from = compartment,
		values_from = Abundance,
		values_fill = list(Abundance = 0)
	) %>% group_by(pathway, species, index) %>%
	summarise(
		bulk_soil = sum(bulk_soil),
		rhizosheath = sum(rhizosheath),
		.groups = 'drop'
	) %>% mutate(
		delta = rhizosheath - bulk_soil,
		FC = ifelse(
			bulk_soil == 0,
			rhizosheath,
			rhizosheath/bulk_soil
		)
	) %>% arrange(desc(FC)) %>%
	mutate( 
		across(
			where(is.numeric),
			~ num(., digits = 6)
		), 
		FC_rank = row_number()
	)

library("file2meco")
library("ALDEx2")

delta.pvals <- delta.paths %>%
	group_by(species, pathway) %>%
	summarise(
		p.val = tryCatch(
			wilcox.test(rhizosheath, bulk_soil, paired = T)$p.value,
			error = function(e) NA_real_
		),
		FC.med = median(as.numeric(rhizosheath) - as.numeric(bulk_soil)),
		.groups = 'drop'
	) %>%
	mutate(
		logFC = ifelse(FC.med == 0, 0, log2(abs(FC.med))),
		logFC = ifelse(FC.med < 0, -logFC, logFC),
		neglog10p = -log10(p.val)
	)

delta.ranx <- delta.paths %>%
	group_by(pathway, species) %>%
	summarise(
		mean_rank = mean(FC_rank, na.rm = T),
		.groups = 'drop'
	) %>% arrange(-desc(mean_rank))

delta.UA <- delta.paths %>%
	group_by(pathway) %>%
	summarise(
		mean_rank = mean(FC_rank, na.rm = T),
		.groups = 'drop'
	) %>% arrange(-desc(mean_rank))

################### VOLCANO-PLOT #######################
# The best input for this, in reality, are the KEGG Orthologues
# Therefore, I'll load those files as input
# Then normalize counts using DESEQ2 of RS vs BS

library("DESeq2")
library("data.table")
KO.counts <- fread("./picrust_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz") %>%
	as.data.frame()
KO.counts <- KO.counts[, sort(names(KO.counts))]
# head(KO.counts)
names(KO.counts)[1] <- "functions"
rownames(KO.counts) <- KO.counts$functions
# head(KO.counts)
KO.counts$functions <- NULL

# all(colnames(KO.counts) == sample.info$Sample)

sample.info$species <- as.factor(sample.info$species)
sample.info$compartment <- as.factor(sample.info$compartment)

new.sp.all.results <- map(levels(sample.info$species), function(sp){
	tmp.sample.info <- sample.info %>% filter(species == sp)
	tmp.KO.counts <- KO.counts[, tmp.sample.info$Sample]
	tmp.KO.counts <- tmp.KO.counts[rowSums(KO.counts[]) > 0, ]
	tmp.deseq <- DESeqDataSetFromMatrix(
		countData = round(tmp.KO.counts*1000),
		colData = tmp.sample.info,
		design = ~ compartment
	)
	tmp.deseq <- estimateSizeFactors(tmp.deseq, type = "poscounts")
	tmp.deseq <- DESeq(tmp.deseq, test = 'Wald', fitType = 'local')
	tmp.results <- results(tmp.deseq,
		contrast = c("compartment", "rhizosheath", "bulk_soil")
	)
	tmp.deseq.df <- as.data.frame(tmp.results) %>%
	rownames_to_column("KOs") %>%
	mutate(
		species = sp,
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
			padj < 0.05 & log2FoldChange > 2 ~ "Up in RS",
			padj < 0.05 & log2FoldChange < -2 ~ "Up in BS",
			TRUE ~ "Non-significant"
		)
	)
	return(tmp.deseq.df)
	}
)

#KO.deseq <- DESeqDataSetFromMatrix(
#		countData = round(KO.counts*1000),
#		colData = sample.info,
#		design = species ~ compartment
#	)
#
#deseq.alt <- estimateSizeFactors(KO.deseq, type = "poscounts")
#
#deseq.alt <- DESeq(deseq.alt, test = 'Wald', fitType = 'local')
#
#deseq.alt.res <-  results(
#		deseq.alt,
#		contrast = c("compartment", "rhizosheath", "bulk_soil")
#	)
#
#deseq.alt.df <- as.data.frame(deseq.alt.res) %>%
#	rownames_to_column("KOs") %>%
#	mutate(
#		log2FoldChange = ifelse(
#			is.na(log2FoldChange),
#			0,
#			log2FoldChange
#		),
#		padj = ifelse(
#			is.na(padj),
#			1,
#			padj
#		),
#		Significance = case_when(
#			padj < 0.05 & log2FoldChange > 20 ~ "Up in RH",
#			padj < 0.05 & log2FoldChange < -20 ~ "Up in BS",
#			TRUE ~ "Non-significant"
#		),
#		Label = ifelse(
#			padj < 0.001 & abs(log2FoldChange) > 20,
#			KOs,
#			NA
#		)
#	)

new.sp.all.results <- bind_rows(new.sp.all.results)
write.csv(sp.all.dfs, "./picrust_out/Volcano_KO_data.csv")

#volcano.alt <- ggplot(
#	deseq.alt.df,
#	aes(
#		x = log2FoldChange, y = -log10(padj),
#		color = Significance,
##		label = Label
#	))+geom_point(alpha = 0.6, size = 2.5)+
#	geom_vline(
#		xintercept = c(-20, 20),
#		linetype = "dashed",
#		alpha = 0.5
#	)+geom_hline(
#		yintercept = -log10(0.05),
#		linetype = "dashed",
#		alpha = 0.5
#	)+scale_color_manual(
#		values = c("Up in RH" = "wheat",
#			"Up in BS" = "brown",
#			"Non-significant" = "gray"
#		)
##	)+geom_text_repel(
##		max.overlaps = 20,
##		size = 3,
##		box.padding = 0.5,
##		segment.color = "grey50"
#	)+labs(
#		x = "log2(Fold Change) Rhizosheath (RH) vs Bulk Soil (BS)",
#		y = "-log10(Adjusted p-value)",
#		title = "Differential Abundance: Compartment RH vs BS",
##		subtitle = "Controlled for species effects (n=4 species)"
#	)+theme_minimal(base_size = 12)+
#	theme(
#		legend.position = "top",
#		panel.grid.minor = element_blank(),
#		plot.title = element_text(face = "bold")
#	)+scale_shape_manual(
#			values = c("High presence in RH" = 24,
#			"Low presence in RH" = 25,
#			"Non-significant" = 21
#		)
#	)
#volcano.alt

library("ggh4x")
volcano.all.sp <- ggplot(
	sp.all.results,
	aes(
		x = log2FoldChange, y = -log10(padj),
		color = Significance
	))+
	geom_point(alpha = 0.6, size = 2.5) +
	geom_vline(xintercept = c(-10, 10), linetype = "dashed", alpha = 0.5) +
	geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
	scale_color_manual(
		values = c("Up in RS" = "wheat",
			"Up in BS" = "brown",
			"Non-significant" = "gray")
	)+
	labs(
		x = "log2(Fold Change) Rhizosheath (RS) vs Bulk Soil (BS)",
		y = "-log10(Adjusted p-value)",
		title = "Differential Abundance: Compartment RS vs BS",
	)+
	theme_minimal(base_size = 12)+
	facet_wrap2(
		~species, scales = 'free_x',
		strip = strip_themed(background_x = elem_list_rect(
			fill = c("#f66c42", "#2995bc", "#43b284", "#fab255")
		))
	)+
	theme(
		legend.position = 'top',
		panel.grid.minor = element_blank(),
		plot.title = element_text(face = "bold", size = 36),
		panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
		strip.text = element_text(face = "bold.italic", color = "black", size = 30),
		strip.background = element_rect(linewidth = 0.8),
		axis.title = element_text(size = 30, face = "bold"),
		axis.text = element_text(size = 24),
		legend.text = element_text(size = 30)
    )

jpeg("./Heatmaps/species_spec_VolcanoPlots.jpeg", height = 4000, width = 4000, res = 300)
volcano.all.sp
dev.off()

new.down.RS <- new.sp.all.results[new.sp.all.results$Significance == "Up in BS", ]
new.list_functions <- split(new.sp.all.results$KOs, new.sp.all.results$species)

require(ggvenn)
new.core.functions <- ggvenn(
		new.list_functions,
		fill_color = c("#f66c42", "#fab255", "#43b284", "#2995bc"),
		text_size = 7,
		set_name_size = 7
	)+
	labs(title = "RS-enriched share functions among MHSC species")+
	theme(plot.title = element_text(
			size = 30,
			face = "bold",
			hjust = 0.5
		)
	)+coord_cartesian(clip = 'off')
core.functions$layers[[3]]$aes_params$fontface <- "bold.italic"
core.functions

core_asvs <- Reduce(intersect, list_omne)
taxa_core <- tax_table(rseq)[core_asvs,]













results(deseq.alt, contrast = c("compartment", "bulk_soil", "rhizosheath")) %>%
	as.data.frame() %>%
	rownames_to_column("KOs") %>%
	mutate(
		log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
		padj = ifelse(is.na(padj), 1, padj),
		Significance = case_when(
			padj < 0.05 & log2FoldChange > 20 ~ "Up in BS",
			padj < 0.05 & log2FoldChange < -20 ~ "Down in BS",
			TRUE ~ "Non-significant"
		),
		Label = ifelse(
			padj < 0.001 & abs(log2FoldChange) > 20,
			KOs,
			NA
		)
	) %>% ggplot(
		aes(x = log2FoldChange, y = -log10(padj), color = Significance))+
	geom_point(alpha = 0.6, size = 2.5)+
	geom_vline(xintercept = c(-20, 20), linetype = "dashed", alpha = 0.5)+
	geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5)+
	scale_color_manual(values = c("Down in BS" = "wheat", "Up in BS" = "brown", "Non-significant" = "gray"))+
	labs(
		x = "log2(Fold Change) Bulk Soil (BS) vs Rhizosheath (RS)",
		y = "-log10(Adjusted p-value)",
		title = "Differential Abundance: Compartment RS vs BS",
	)+theme_minimal(base_size = 12)+
	theme(
		legend.position = "top",
		panel.grid.minor = element_blank(),
		plot.title = element_text(face = "bold")
	)+scale_shape_manual(values = c("Up in BS" = 25, "Down in BS" = 24, "Non-significant" = 21
		)
	)
volcano.alt



volcano.deseq <- ggplot(
		deseq.res,
		aes(
			x = log2FoldChange,
			y = -log10(padj)
		)
	)+
	geom_point()+
	labs(
		x = "Log2 Fold Change",
		y = "-Log10 adjusted p-value"
	)

volcano.deseq


volcano.overall <- delta.pvals %>%
	mutate(
		log2FC = log2(FC + 1e-8),
		neglog10delta = -log10(abs(delta) + 1)
	)

ggplot(delta.pvals,
	aes(
		x = logFC,
		y = neglog10p
	)
)+
geom_point(size = 1.2, col = "grey20")

library("file2meco")


library(ggh4x)
jpeg("./Heatmaps/barplot_top_pathways.jpeg", height = 4000, width = 4000, res = 300)
top.plot <- ggplot(top.paths,
		aes(x = mean_abund, y = pathway, fill = compartment)
	)+geom_col(position = position_dodge(width = 0.8))+
	geom_errorbar(
		aes(
			xmin = mean_abund - sd_abund,
			xmax = mean_abund + sd_abund
		),
		width = 0.3,
		position = position_dodge(width = 0.8)
	)+facet_wrap2(
		~species,
		scales = 'free_y',
		strip = strip_themed(
			background_x = elem_list_rect(fill = sp.cols)
		), labeller = label_wrap_gen(width = 13)
	)+
	labs(
		x = "Predicted abundance",
		y = NULL,
		fill = "Compartment"
	)+
	theme_bw(base_size = 14)+
	theme(strip.text.x = element_text(size = 12, face = 'bold.italic'))+
	scale_fill_manual(
		values = c(
			"bulk_soil" = 'sienna',
			"rhizosheath" = 'wheat'
		), labels = c(
			"bulk_soil" = "Bulk Soil",
			"rhizosheath" = "Rhizosheath"
		)
	)
top.plot
dev.off()

################### VOLCANO-PLOT ATTEMPT 2 #######################
# The best input for this, in reality, are the KEGG Orthologues
# Therefore, I'll load those files as input
# Then normalize counts using DESEQ2 of RS vs BS

library("DESeq2")
library("data.table")
KO.counts <- fread("./picrust_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz") %>%
	as.data.frame()
KO.counts <- KO.counts[, sort(names(KO.counts))]
# head(KO.counts)
names(KO.counts)[1] <- "functions"
rownames(KO.counts) <- KO.counts$functions
# head(KO.counts)
KO.counts$functions <- NULL

# all(colnames(KO.counts) == sample.info$Sample)

sample.info$species <- as.factor(sample.info$species)
sample.info$compartment <- as.factor(sample.info$compartment)
new.sample.info <- sample.info
new.sample.info$group <- factor(paste0(sample.info$species,"_",sample.info$compartment))

global.deseq <- DESeqDataSetFromMatrix(
		countData = round(KO.counts*100),
		colData = new.sample.info,
		design = ~ group
	)
global.deseq <- estimateSizeFactors(global.deseq, type = "poscounts")
global.deseq <- DESeq(global.deseq, test = 'Wald', fitType = 'local')

new.sp.all.results <- map(levels(new.sample.info$species), function(sp){
	tmp.results <- results(global.deseq,
		contrast = c("group", paste0(sp,"_rhizosheath"), paste0(sp,"_bulk_soil"))
	)
	tmp.deseq.df <- as.data.frame(tmp.results) %>%
	rownames_to_column("KOs") %>%
	mutate(
		species = sp,
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
			padj < 0.05 & log2FoldChange > 2 ~ "Up in RS",
			padj < 0.05 & log2FoldChange < -2 ~ "Up in BS",
			TRUE ~ "Non-significant"
		)
	)
	return(tmp.deseq.df)
	}
)
str(new.sp.all.results)

new.sp.all.results <- bind_rows(new.sp.all.results)
write.csv(new.sp.all.results, "./picrust_out/narB_Volcano_KO_data.csv")

require("ggh4x")
narB.volcano.all.sp <- ggplot(
	new.sp.all.results,
	aes(
		x = log2FoldChange, y = -log10(padj),
		color = Significance
	))+
	geom_point(alpha = 0.6, size = 2.5) +
	geom_vline(xintercept = c(-2, 2), linetype = "dashed", alpha = 0.5) +
	geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
	scale_color_manual(
		values = c("Up in RS" = "wheat",
			"Up in BS" = "brown",
			"Non-significant" = "gray")
	)+
	labs(
		x = "log2(Fold Change) Rhizosheath (RS) vs Bulk Soil (BS)",
		y = "-log10(Adjusted p-value)",
		title = "Differential Abundance: Compartment RS vs BS",
	)+
	theme_minimal(base_size = 12)+
	facet_wrap2(
		~species,
		strip = strip_themed(background_x = elem_list_rect(
			fill = c("#f66c42", "#2995bc", "#43b284", "#fab255")
		))
	)+
	theme(
		legend.position = 'top',
		panel.grid.minor = element_blank(),
		plot.title = element_text(face = "bold", size = 24),
		panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
		strip.text = element_text(face = "bold.italic", color = "black", size = 18),
		strip.background = element_rect(linewidth = 0.8),
		axis.title = element_text(size = 18, face = "bold"),
		axis.text = element_text(size = 15),
		legend.text = element_text(size = 15)
    )

jpeg("./Heatmaps/narB_species_spec_VolcanoPlots.jpeg", height = 4000, width = 4000, res = 300)
narB.volcano.all.sp
dev.off()

new.up.RS <- new.sp.all.results[new.sp.all.results$Significance == "Up in RS", ]
new.list_functions <- split(new.up.RS$KOs, new.up.RS$species)

require(ggvenn)
narB.core.functions <- ggvenn(
		new.list_functions,
		fill_color = c("#2995bc", "#fab255", "#f66c42", "#43b284"),
		text_size = 7,
		set_name_size = 7
	)+
	labs(title = "RS-enriched share functions among MHSC species")+
	theme(plot.title = element_text(
			size = 30,
			face = "bold",
			hjust = 0.5
		)
	)+coord_cartesian(clip = 'off')
narB.core.functions$layers[[3]]$aes_params$fontface <- "italic"
narB.core.functions

narB.list.functions <- Reduce(intersect, new.list_functions)
taxa_core <- tax_table(rseq)[core_asvs,]













results(deseq.alt, contrast = c("compartment", "bulk_soil", "rhizosheath")) %>%
	as.data.frame() %>%
	rownames_to_column("KOs") %>%
	mutate(
		log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
		padj = ifelse(is.na(padj), 1, padj),
		Significance = case_when(
			padj < 0.05 & log2FoldChange > 20 ~ "Up in BS",
			padj < 0.05 & log2FoldChange < -20 ~ "Down in BS",
			TRUE ~ "Non-significant"
		),
		Label = ifelse(
			padj < 0.001 & abs(log2FoldChange) > 20,
			KOs,
			NA
		)
	) %>% ggplot(
		aes(x = log2FoldChange, y = -log10(padj), color = Significance))+
	geom_point(alpha = 0.6, size = 2.5)+
	geom_vline(xintercept = c(-20, 20), linetype = "dashed", alpha = 0.5)+
	geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5)+
	scale_color_manual(values = c("Down in BS" = "wheat", "Up in BS" = "brown", "Non-significant" = "gray"))+
	labs(
		x = "log2(Fold Change) Bulk Soil (BS) vs Rhizosheath (RS)",
		y = "-log10(Adjusted p-value)",
		title = "Differential Abundance: Compartment RS vs BS",
	)+theme_minimal(base_size = 12)+
	theme(
		legend.position = "top",
		panel.grid.minor = element_blank(),
		plot.title = element_text(face = "bold")
	)+scale_shape_manual(values = c("Up in BS" = 25, "Down in BS" = 24, "Non-significant" = 21
		)
	)
volcano.alt



volcano.deseq <- ggplot(
		deseq.res,
		aes(
			x = log2FoldChange,
			y = -log10(padj)
		)
	)+
	geom_point()+
	labs(
		x = "Log2 Fold Change",
		y = "-Log10 adjusted p-value"
	)

volcano.deseq


volcano.overall <- delta.pvals %>%
	mutate(
		log2FC = log2(FC + 1e-8),
		neglog10delta = -log10(abs(delta) + 1)
	)

ggplot(delta.pvals,
	aes(
		x = logFC,
		y = neglog10p
	)
)+
geom_point(size = 1.2, col = "grey20")

library("file2meco")


library(ggh4x)
jpeg("./Heatmaps/barplot_top_pathways.jpeg", height = 4000, width = 4000, res = 300)
top.plot <- ggplot(top.paths,
		aes(x = mean_abund, y = pathway, fill = compartment)
	)+geom_col(position = position_dodge(width = 0.8))+
	geom_errorbar(
		aes(
			xmin = mean_abund - sd_abund,
			xmax = mean_abund + sd_abund
		),
		width = 0.3,
		position = position_dodge(width = 0.8)
	)+facet_wrap2(
		~species,
		scales = 'free_y',
		strip = strip_themed(
			background_x = elem_list_rect(fill = sp.cols)
		), labeller = label_wrap_gen(width = 13)
	)+
	labs(
		x = "Predicted abundance",
		y = NULL,
		fill = "Compartment"
	)+
	theme_bw(base_size = 14)+
	theme(strip.text.x = element_text(size = 12, face = 'bold.italic'))+
	scale_fill_manual(
		values = c(
			"bulk_soil" = 'sienna',
			"rhizosheath" = 'wheat'
		), labels = c(
			"bulk_soil" = "Bulk Soil",
			"rhizosheath" = "Rhizosheath"
		)
	)
top.plot
dev.off()

################### VOLCANO-PLOT ATTEMPT 3 #######################
# The best input for this, in reality, are the KEGG Orthologues
# Therefore, I'll load those files as input
# Then normalize counts using ALDEx2 of RS vs BS (after Nearing et. al, 2022)

library("ALDEx2")
library("data.table")
KO.counts <- fread("./picrust_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz") %>%
	as.data.frame()
KO.counts <- KO.counts[, sort(names(KO.counts))]
# head(KO.counts)
names(KO.counts)[1] <- "functions"
rownames(KO.counts) <- KO.counts$functions
# head(KO.counts)
KO.counts$functions <- NULL

# all(colnames(KO.counts) == sample.info$Sample)

sample.info$species <- as.factor(sample.info$species)
sample.info$compartment <- as.factor(sample.info$compartment)
Ax2.sample.info <- sample.info
Ax2.sample.info$group <- factor(paste0(sample.info$species,"_",sample.info$compartment))

global.aldex <- aldex.clr(
		reads = round(as.matrix(KO.counts)*100),
		conds = as.character(Ax2.sample.info$compartment),
		mc.sample = 150,
		denom = "all",
		verbose = T
	)

save.image("./RSessions/Jan-12.RData")

aldex.sp.all.results <- map(
	levels(Ax2.sample.info$species),
	function(sp){
		tmp.id <- Ax2.sample.info$species == sp
		tmp.counts <- KO.counts[, tmp.id]
		tmp.conds <- as.character(Ax2.sample.info$compartment[tmp.id])
		tmp.aldex <- aldex.clr(
			reads = round(as.matrix(tmp.counts)*100),
			conds = tmp.conds,
			mc.sample = 128,
			denom = 'all',
			verbose = T
		)
		tmp.tt <- aldex.ttest(tmp.aldex)
		tmp.ef <- aldex.effect(tmp.aldex)
		aldex.res <- cbind(tmp.tt, tmp.ef) %>%
			as.data.frame() %>%
			rownames_to_column("KOs") %>%
			mutate(
				species = sp,
				effect = ifelse(is.na(effect), 0, effect), #Log2-Fold Change 	
				we.ep = ifelse(is.na(we.ep), 1, we.ep), #p-value
				Significance = case_when(
					we.ep < 0.05 & effect > 1 ~ "Up in RS",
					we.ep < 0.05 & effect < -1 ~ "Up in BS",
					TRUE ~ "Non-significant"
				)
			)
		return(aldex.res)
		}
	)
str(aldex.sp.all.results)

aldex.sp.all.results <- bind_rows(aldex.sp.all.results)
write.csv(aldex.sp.all.results, "./picrust_out/ALDEx2_Volcano_KO_data.csv")

require("ggh4x")
aldex.volcano.all.sp__old <- aldex.volcano.all.sp
aldex.volcano.all.sp <- ggplot(
	aldex.sp.all.results,
	aes(
		x = effect, y = -log10(we.ep),
		color = Significance
	))+
	geom_point(alpha = 0.6, size = 2.5) +
	geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
	geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
	scale_color_manual(
		values = c("Up in RS" = "wheat",
			"Up in BS" = "brown",
			"Non-significant" = "gray")
	)+
	labs(
		x = "ALDEx2 effect size (RS vs BS)",
		y = "-log10 (BH-adjusted expected p-value)",
		title = "Differential Abundance (ALDEx2): Rhizosheath vs Bulk Soil",
	)+
	theme_minimal(base_size = 12)+
	facet_wrap2(
		~species,
		strip = strip_themed(background_x = elem_list_rect(
			fill = c("#f66c42", "#2995bc", "#43b284", "#fab255")
		))
	)+
	theme(
		legend.position = 'top',
		panel.grid.minor = element_blank(),
		plot.title = element_text(hjust = 0.5, face = 'bold', size = 32),
		panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.8),
		strip.text = element_text(face = 'bold.italic', color = 'black', size = 26),
		strip.background = element_rect(linewidth = 0.8),
		axis.title = element_text(size = 26, face = 'bold'),
		axis.text = element_text(size = 20),
		legend.text = element_text(size = 26),
		legend.title = element_text(size = 26, face = 'bold')
    )

jpeg("./Heatmaps/ALDEx2_species_VolcanoPlots_v2.jpeg", height = 8000, width = 4000, res = 300)
aldex.volcano.all.sp
dev.off()

aldex.up.RS <- aldex.sp.all.results[aldex.sp.all.results$Significance == "Up in RS", ]
aldex.list_functions <- split(aldex.up.RS$KOs, aldex.up.RS$species)

require(ggvenn)
aldex.core.functions <- ggvenn(
		aldex.list_functions,
		fill_color = c("#f66c42", "#2995bc", "#43b284", "#fab255"),
		text_size = 7,
		set_name_size = 7
#		set_name_color = "white"
	)+
	labs(title = "RS-enriched share functions among MHSC species")+
	theme(
		plot.title = element_text(size = 30, face = "bold", hjust = 0.5)
#		panel.background = element_rect(fill = "black")
	)+coord_cartesian(clip = 'off')
aldex.core.functions$layers[[3]]$aes_params$fontface <- "bold.italic"

jpeg("./Heatmaps/ALDEx2_species_FunctionalCore.jpeg", height = 4000, width = 4000, res = 300)
aldex.core.functions
dev.off()

aldex.list.core_functions <- Reduce(intersect, aldex.list_functions)

################### VOLCANO-PLOT OF TAXA_CORE'S KOs #######################
library("ALDEx2")
library("data.table")
KO.strats <- fread("./picrust_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz") %>%
	as.data.frame()
rownames(KO.strats) <- KO.strats$sequence
KO.strats <- KO.strats[,-1]
# all(core_asvs %in% rownames(KO.strats))
KO.strats.core <- KO.strats[intersect(core_asvs, rownames(KO.strats)), ]
KO.strats.core <- KO.strats.core[, colSums(KO.strats.core) > 0]

filt.positions <- which(colnames(KO.strats.core) %in% aldex.up.RS$KOs)

KO.strats.filt <- KO.strats.core[,filt.positions]
KO.strats.filt_core <- KO.strats.filt[,which(colnames(KO.strats.filt) %in% aldex.list.core_functions)]

# aldex.list.core_functions[-2] %in% names(KO.strats.filt_core) #Sanity check
KO.strats.filt_core <- cbind(KO.strats.filt_core, "ko:K00846" = rep(0, 17))

#Creating a metadata tracking dataframe for KO terms
abbs_vector <- c("ko:K00842", "KHK", "atoD", "ybtX", "cidA", "iolT", "ydbD",
	"araJ", "ko:K08989", "psuK", "adrA", "bar", "ladA", "fosX")
functions_vector <- c("Unknown function", "Plant-Fungi signalling",
	"Exudated C recyclying", "Fe and Zn scavenging", "Stress responses",
	"Exudated C recyclying", "Stress responses", "Exudated C recyclying",
	"Unknown function", "Necromass N cyclying", "Stress responses",
	"Necromass N cyclying", "Necromass N cyclying", "Stress responses")
visual_mapping_labels <- as.data.frame(
	cbind("KO_terms" = sort(names(KO.strats.filt_core)),
		"Abbreviations" = abbs_vector,
		"Functional_PWs" = functions_vector)
	) |> arrange(tolower(Functional_PWs))
for(i in 1:14){
	if(visual_mapping_labels[i,1] == visual_mapping_labels[i,2]){
		visual_mapping_labels[i, 4] <- visual_mapping_labels[i,2]
	} else {
		visual_mapping_labels[i, 4] <- paste0(
			visual_mapping_labels[i,2],
			" (",
			visual_mapping_labels[i,1],
			")"
		)
	}
}
names(visual_mapping_labels)[4] <- "Notation"
# Creating the heatmap
BiocManager::install("ComplexHeatmap")
require("ComplexHeatmap")
#rows_hclust <- hclust(dist(KO.strats.filt_bin)) # Create Clustering for rows
#order.rows_hclust <- reorder(as.dendrogram(rows_hclust), rowSums(KO.strats.filt_bin))
cols_FPW <- setNames(c("purple", "indianred", "aquamarine", "magenta", "darkred", "gray60"),
	unique(visual_mapping_labels$Functional_PWs))
col_not <- HeatmapAnnotation( # Map colours according to function
	Pathway = visual_mapping_labels$Functional_PWs,
	col = list(Pathway = cols_FPW),
	annotation_legend_param = list(Pathway = 
		list(
			nrow = 2,
			title_gp = gpar(fontsize = 20, fontface = 'bold'),
			labels_gp = gpar(fontsize = 18)
		)
	),
	show_annotation_name = FALSE
)
# Map ASV and Tax labels to matrix
Tax_ASV_lab <- prune.table |> as.data.frame() |> rownames_to_column() |> 
	separate(col = rowname, into = c("Tax_lab", "ASV_lab"), sep = ' : ') |> 
	select(Tax_lab, ASV_lab) |> arrange(ASV_lab) |> 
	unite("Tax_ASV_lab", Tax_lab, ASV_lab, sep = ' : ')
# Create binary matrix
KO.strats.filt_bin <- as.data.frame(ifelse(KO.strats.filt_core > 0, 1, 0)) |> 
	rownames_to_column() |> arrange(rowname) |> mutate(Tax_ASV_lab) |> 
	column_to_rownames("Tax_ASV_lab") |> select(-rowname) |> as.matrix()
KO.strats.filt_bin <- KO.strats.filt_bin[,visual_mapping_labels$KO_terms]
jpeg("./Heatmaps/Intersection_Cores_v2.jpeg", height = 4000, width = 4000, res = 300)
draw(
	Heatmap(
		KO.strats.filt_bin,
		name = NULL,
		show_heatmap_legend = FALSE,
		col = c("white", "green"),
		cluster_columns = FALSE,
		cluster_rows = TRUE,
		row_labels = gsub(' : ', '\n', rownames(KO.strats.filt_bin)),
		column_labels = visual_mapping_labels$Notation,
		top_annotation = col_not,
		row_names_gp = gpar(fontsize = 16),
		column_names_gp = gpar(fontsize = 18, fontface = 'bold'),
		rect_gp = gpar(type = 'none'),
		cell_fun = function(j, i, x, y, width, height, fill){
			if(KO.strats.filt_bin[i, j] == 1){
				grid.circle(x = x, y = y,
					r = min(unit.c(width, height))*0.5,
					gp = gpar(fill = 'green', col = NA)
				)
			} else {
				grid.circle(x = x, y = y,
					r = min(unit.c(width, height))*0.25,
					gp = gpar(fill = 'gray90', col = NA)
				)
			}
		}
	), annotation_legend_side = "top",
	padding = unit(c(2, 2, 2, 57), 'mm')
)
dev.off()

sample.info$species <- as.factor(sample.info$species)
sample.info$compartment <- as.factor(sample.info$compartment)
Ax2.sample.info <- sample.info
Ax2.sample.info$group <- factor(paste0(sample.info$species,"_",sample.info$compartment))

global.aldex <- aldex.clr(
		reads = round(as.matrix(KO.counts)*100),
		conds = as.character(Ax2.sample.info$compartment),
		mc.sample = 150,
		denom = "all",
		verbose = T
	)

save.image("./RSessions/Jan-12.RData")

aldex.sp.all.results <- map(
	levels(Ax2.sample.info$species),
	function(sp){
		tmp.id <- Ax2.sample.info$species == sp
		tmp.counts <- KO.counts[, tmp.id]
		tmp.conds <- as.character(Ax2.sample.info$compartment[tmp.id])
		tmp.aldex <- aldex.clr(
			reads = round(as.matrix(tmp.counts)*100),
			conds = tmp.conds,
			mc.sample = 128,
			denom = 'all',
			verbose = T
		)
		tmp.tt <- aldex.ttest(tmp.aldex)
		tmp.ef <- aldex.effect(tmp.aldex)
		aldex.res <- cbind(tmp.tt, tmp.ef) %>%
			as.data.frame() %>%
			rownames_to_column("KOs") %>%
			mutate(
				species = sp,
				effect = ifelse(is.na(effect), 0, effect), #Log2-Fold Change 	
				we.ep = ifelse(is.na(we.ep), 1, we.ep), #p-value
				Significance = case_when(
					we.ep < 0.05 & effect > 1 ~ "Up in RS",
					we.ep < 0.05 & effect < -1 ~ "Up in BS",
					TRUE ~ "Non-significant"
				)
			)
		return(aldex.res)
		}
	)
str(aldex.sp.all.results)

aldex.sp.all.results <- bind_rows(aldex.sp.all.results)
write.csv(aldex.sp.all.results, "./picrust_out/ALDEx2_Volcano_KO_data.csv")

require("ggh4x")
aldex.volcano.all.sp <- ggplot(
	aldex.sp.all.results,
	aes(
		x = effect, y = -log10(we.ep),
		color = Significance
	))+
	geom_point(alpha = 0.6, size = 2.5) +
	geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
	geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
	scale_color_manual(
		values = c("Up in RS" = "wheat",
			"Up in BS" = "brown",
			"Non-significant" = "gray")
	)+
	labs(
		x = "ALDEx2 effect size (RS vs BS)",
		y = "-log10 (BH-adjusted expected p-value)",
		title = "Differential Abundance (ALDEx2): Rhizosheath vs Bulk Soil",
	)+
	theme_minimal(base_size = 12)+
	facet_wrap2(
		~species,
		strip = strip_themed(background_x = elem_list_rect(
			fill = c("#f66c42", "#2995bc", "#43b284", "#fab255")
		))
	)+
	theme(
		legend.position = 'top',
		panel.grid.minor = element_blank(),
		plot.title = element_text(face = "bold", size = 24),
		panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
		strip.text = element_text(face = "bold.italic", color = "black", size = 18),
		strip.background = element_rect(linewidth = 0.8),
		axis.title = element_text(size = 18, face = "bold"),
		axis.text = element_text(size = 15),
		legend.text = element_text(size = 15)
    )

jpeg("./Heatmaps/ALDEx2_species_VolcanoPlots.jpeg", height = 4000, width = 4000, res = 300)
aldex.volcano.all.sp
dev.off()

aldex.up.RS <- aldex.sp.all.results[aldex.sp.all.results$Significance == "Up in RS", ]
aldex.list_functions <- split(aldex.up.RS$KOs, aldex.up.RS$species)

require(ggvenn)
aldex.core.functions <- ggvenn(
		aldex.list_functions,
		fill_color = c("#2995bc", "#fab255", "#f66c42", "#43b284"),
		text_size = 7,
		set_name_size = 7
	)+
	labs(title = "RS-enriched share functions among MHSC species")+
	theme(plot.title = element_text(
			size = 30,
			face = "bold",
			hjust = 0.5
		)
	)+coord_cartesian(clip = 'off')
aldex.core.functions$layers[[3]]$aes_params$fontface <- "italic"
aldex.core.functions

aldex.list.core_functions <- Reduce(intersect, aldex.list_functions)
taxa_core <- tax_table(rseq)[core_asvs,]















results(deseq.alt, contrast = c("compartment", "bulk_soil", "rhizosheath")) %>%
	as.data.frame() %>%
	rownames_to_column("KOs") %>%
	mutate(
		log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
		padj = ifelse(is.na(padj), 1, padj),
		Significance = case_when(
			padj < 0.05 & log2FoldChange > 20 ~ "Up in BS",
			padj < 0.05 & log2FoldChange < -20 ~ "Down in BS",
			TRUE ~ "Non-significant"
		),
		Label = ifelse(
			padj < 0.001 & abs(log2FoldChange) > 20,
			KOs,
			NA
		)
	) %>% ggplot(
		aes(x = log2FoldChange, y = -log10(padj), color = Significance))+
	geom_point(alpha = 0.6, size = 2.5)+
	geom_vline(xintercept = c(-20, 20), linetype = "dashed", alpha = 0.5)+
	geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5)+
	scale_color_manual(values = c("Down in BS" = "wheat", "Up in BS" = "brown", "Non-significant" = "gray"))+
	labs(
		x = "log2(Fold Change) Bulk Soil (BS) vs Rhizosheath (RS)",
		y = "-log10(Adjusted p-value)",
		title = "Differential Abundance: Compartment RS vs BS",
	)+theme_minimal(base_size = 12)+
	theme(
		legend.position = "top",
		panel.grid.minor = element_blank(),
		plot.title = element_text(face = "bold")
	)+scale_shape_manual(values = c("Up in BS" = 25, "Down in BS" = 24, "Non-significant" = 21
		)
	)
volcano.alt



volcano.deseq <- ggplot(
		deseq.res,
		aes(
			x = log2FoldChange,
			y = -log10(padj)
		)
	)+
	geom_point()+
	labs(
		x = "Log2 Fold Change",
		y = "-Log10 adjusted p-value"
	)

volcano.deseq


volcano.overall <- delta.pvals %>%
	mutate(
		log2FC = log2(FC + 1e-8),
		neglog10delta = -log10(abs(delta) + 1)
	)

ggplot(delta.pvals,
	aes(
		x = logFC,
		y = neglog10p
	)
)+
geom_point(size = 1.2, col = "grey20")

library("file2meco")


library(ggh4x)
jpeg("./Heatmaps/barplot_top_pathways.jpeg", height = 4000, width = 4000, res = 300)
top.plot <- ggplot(top.paths,
		aes(x = mean_abund, y = pathway, fill = compartment)
	)+geom_col(position = position_dodge(width = 0.8))+
	geom_errorbar(
		aes(
			xmin = mean_abund - sd_abund,
			xmax = mean_abund + sd_abund
		),
		width = 0.3,
		position = position_dodge(width = 0.8)
	)+facet_wrap2(
		~species,
		scales = 'free_y',
		strip = strip_themed(
			background_x = elem_list_rect(fill = sp.cols)
		), labeller = label_wrap_gen(width = 13)
	)+
	labs(
		x = "Predicted abundance",
		y = NULL,
		fill = "Compartment"
	)+
	theme_bw(base_size = 14)+
	theme(strip.text.x = element_text(size = 12, face = 'bold.italic'))+
	scale_fill_manual(
		values = c(
			"bulk_soil" = 'sienna',
			"rhizosheath" = 'wheat'
		), labels = c(
			"bulk_soil" = "Bulk Soil",
			"rhizosheath" = "Rhizosheath"
		)
	)
top.plot
dev.off()