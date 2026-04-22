##### ESTE SCRIPT INDICA EL PROCESO PARA ELABORAR LOS DIAGRAMAS DE VENN ####
require("microbiome")
require("microbiomeutilities")
require("phyloseq")
require("tidyverse")
require("vegan")
require("dplyr")

#### REQUISITOS ####
# Un objeto `phyloseq` con mínimo una tabla de taxonomía y una de abundancias

### Ocurrencia: Mínimo del 0.1% en todas las muestras
rseq <- transform_sample_counts(pseq, function(x) x/sum(x))

#Colapsar columnas de acuerdo al número de ASVs
elbow_ps <- function(ps, rank){ #Function to agglomerate by desired taxonomic rank
	res.bow <- ps %>%
		tax_glom(taxrank = rank) %>% # agglomerate at phylum level
		psmelt() # Melt to long format
	res.bow <- select(res.bow, OTU, Sample, Abundance, species, compartment,
			c(which(colnames(res.bow) == "Phylum"):ncol(res.bow))
		)
	res.bow <- res.bow %>%
		group_by(species, compartment, {{rank}}) %>%
		summarise(Abundance = sum(Abundance), .groups = "drop") %>%
		group_by(species, compartment) %>%
		mutate(Abundance = Abundance/sum(Abundance)) %>%
		ungroup()
	return(res.bow)
}
order.bow <- elbow_ps(pseq, "Order")

host.trts <- unique(meta(rseq)$species)
comp.trts <- unique(meta(rseq)$compartment)
list_omne <- c()
for(n in host.trts){
	m.sub <- subset_samples(rseq, species == n)
	m.sub <- subset_samples(m.sub, compartment == comp.trts[2]) #MR
	m.sub <- prune_taxa(apply(
		otu_table(m.sub),
		1,
		function(x) any(x > 0.001)),
	m.sub)
	omne.m <- rownames(otu_table(m.sub))
#	print(length(rownames(otu_table(m.sub))))
	filt.cut_off <- rowSums(m.sub@otu_table > 0) >= 0.8*nrow(m.sub@sam_data)
#	print(sum(filt.cut_off))
	omne.m <- omne.m[filt.cut_off]
	name <- paste0(n)
	print(paste0("No. of core taxa in ", name, " Rhizosheath: ", length(omne.m)))
	list_omne[[name]] <- omne.m
}

sp.cols <- MetBrewer::met.brewer(name = "Egypt", n = 4, type = "discrete")
require(ggvenn)
core.list <- ggvenn(
		list_omne,
		fill_color = c("#2995bc", "#fab255", "#f66c42", "#43b284"),
		text_size = 7,
		set_name_size = 7
	)+
	labs(title = "RS shared core taxa among MHSC species")+
	theme(plot.title = element_text(
		size = 30,
		face = "bold",
		hjust = 0.5
		),
	)+coord_cartesian(clip = 'off')
core.list$layers[[3]]$aes_params$fontface <- "italic"
core.list

core_asvs <- Reduce(intersect, list_omne)
taxa_core <- tax_table(rseq)[core_asvs,]

per_core_size <- function(pseq, n_perm = 999){
	meta_data <- pseq@sam_data
	observed <- length(core_asvs)
	null.dist <- numeric(n_perm)
	for(i in 1:n_perm){
		shuffled <- sample(meta_data$species)
		tmp.meta <- meta_data
		tmp.meta$species <- shuffled
		tmp.pseq <- pseq
		sample_data(tmp.pseq) <- sample_data(tmp.meta)
		list_null <- list()
		for(n in unique(shuffled)){
			null.sub <- subset_samples(
				tmp.pseq,
				species == n & compartment == "rhizosphere"
			)
			null.sub <- prune_taxa(
				taxa_sums(null.sub) > 0,
				null.sub
			)
			if(ntaxa(null.sub) > 0){
				filt.abd <- apply(
					otu_table(null.sub),
					1,
					function(x) any(x > 0.0001)
				)
				null.sub <- prune_taxa(filt.abd, null.sub)
				filt.occ <- rowSums(otu_table(null.sub) > 0) >= ceiling(0.8*5)
				tmp.core.taxa <- rownames(otu_table(null.sub))[filt.occ]
				list_null[[n]] <- tmp.core.taxa
			}
		}
		if(length(list_null) > 0){
			shared.null <- Reduce(intersect, list_null)
			null.dist[i] <- length(shared.null)
		} else {
			null.dist[i] <- 0
		}
	}
	p_value <- sum(null.dist >= observed + 1) / (n_perm + 1)
	return(list(
		observed = observed,
		p_value = p_value,
		null_dist = null.dist
	))
}

set.seed(120)
p.4999.size <- per_core_size(pseq, n_perm = 4999)

per_core_tags <- function(pseq, n_perm = 999){
	meta.data <- pseq@sam_data
	observed <- core_asvs
	null.dist <- rep(NA, n_perm)
	asv.freq <- setNames(rep(0, length(observed)), observed)
	for(i in 1:n_perm){
		shuffled <- sample(meta.data$species)
		tmp.meta <- meta.data
		tmp.meta$species <- shuffled
		tmp.pseq <- pseq
		sample_data(tmp.pseq) <- sample_data(tmp.meta)
		list_null <- list()
		valid.sp <- 0
		for(n in unique(shuffled)){
			null.sub <- tryCatch({
				subset_samples(
					tmp.pseq,
					species == n & compartment == "rhizosphere"
				)
			}, error = function(e) NULL)
			if(!is.null(null.sub) && nsamples(null.sub) > 0){
				null.sub <- prune_taxa(
					taxa_sums(null.sub) > 0,
					null.sub
				)
				if(ntaxa(null.sub) > 0){
					filt.abd <- apply(
						otu_table(null.sub),
						1,
						function(x) any(x > 0.0001)
					)
					null.sub <- prune_taxa(filt.abd, null.sub)
					if(ntaxa(null.sub) > 0){
						filt.occ <- rowSums(otu_table(null.sub) > 0) >= ceiling(0.8*5)
						if(sum(filt.occ) > 0){
							tmp.core.taxa <- rownames(otu_table(null.sub))[filt.occ]
							list_null[[n]] <- tmp.core.taxa
							valid.sp <- valid.sp + 1
						}
					}
				}
			}
		}
		if(valid.sp >= 2){
			permuted <- tryCatch({
				Reduce(intersect, list_null)
			}, error = function(e) character(0))
#			null.dist[i] <- length(permuted)
			for(asv in observed){
				if(asv %in% permuted){
					asv.freq[asv] <- asv.freq[asv]+1
				}
			}
			tmp.inter <- length(intersect(observed, permuted))
			union.size <- length(union(observed, permuted))
			if(union.size > 0){
				null.dist[i] <- tmp.inter/union.size
			} else {
				null.dist[i] <- 0
			}
		} else {
			null.dist[i] <- NA
		}
	}
	valid.null <- null.dist[!is.na(null.dist)]
	n.valid <- length(valid.null)
	asv.p.val <- sapply(asv.freq, function(freq){
		(sum(freq <= (asv.freq[names(freq)])) + 1) / (n.valid + 1)
	})
	p.val <- (sum(valid.null >= 1) + 1) / (n.valid + 1)
	return(list(
		asv_frequency = asv.freq,
		asv_p_values = asv.p.val,
		null_dist = null.dist,
		general_p_value = p.val,
		valid_permutations = n.valid
	))
}

p.4999.tags <- per_core_tags(pseq = pseq, n = 999)