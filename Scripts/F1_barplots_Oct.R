library("ggtext")

df.order <- psmelt(glom.order)
df.order$Order <- as.factor(df.order$Order)
levels.order <- c(
	levels(df.order$Order)[levels(df.order$Order) != "Other"],
	"Other"
)
df.order$Order <- factor(
	df.order$Order,
	levels = c(
		levels(df.order$Order)[levels(df.order$Order) != "Other"],
		"Other"
	)
)

bar.order <- ggplot(df.order, aes(x = Sample, y = Abundance, fill = Order))+
geom_bar(
	aes(fill = Order),
	stat = 'identity',
	position = "stack"
)+
labs(y = "Relative frequency", x = element_blank())+
theme_minimal()+
theme(
	axis.title = element_text(size = 20),
	axis.text.x = element_blank(),
	legend.key.size = unit(3, 'mm'),
	legend.title = element_text(size = 16, face = 'bold')
)+facet_wrap2(
	~species*compartment,
	ncol = 8,
	scales = 'free',
	strip = strip_themed(
		background_x = elem_list_rect(
			fill = rep(sp.cols, each = 2)
		) 
	), labeller = labeller(
		species = function(x) {
			gsub("(.+)", "<i><b>\\1</b></i>", x) # species in bold italic
		}, compartment = function(x) {
			x <- gsub("bulk_soil", "BS", x)
			x <- gsub("rhizosheath", "RS", x)
			gsub("(.+)", "<b>\\1</b>", x)  # compartment in bold
		}
	)
)+theme(
	axis.text.y = element_blank(),
	axis.ticks.y = element_blank(),
	legend.key.width = unit(0.5, 'cm'),
	legend.text = element_text(size = 10),
	strip.text = element_markdown(
		size = 13.25,
		face = 'bold.italic',
		margin = margin(t = 10),
		lineheight = 1.1
	) 
)+scale_fill_manual(values = randomcoloR::distinctColorPalette(
	k = nrow(NU.order$summary)+1
	)
)+guides(fill = guide_legend(ncol = 1))
bar.order

df.family <- psmelt(glom.family)
df.family$Family <- as.factor(df.family$Family)
df.family$Family <- factor(
	df.family$Family,
	levels = c(
		levels(df.family$Family)[levels(df.family$Family) != "Other"],
		"Other"
	)
)

bar.family <- ggplot(df.family, aes(x = Sample, y = Abundance, fill = Family))+
geom_bar(aes(fill = Family),
	stat = 'identity',
	position = "stack"
)+
labs(y = "Relative frequency", x = element_blank())+
theme_minimal()+
theme(
	axis.title = element_text(size = 20),
	axis.text.x = element_text(size = 11, face = 'bold', angle = 90, hjust = 0),,
	legend.key.size = unit(3, 'mm'),
	legend.title = element_text(size = 16, face = 'bold')
)+facet_wrap2(
	~species*compartment,
	ncol = 8,
	scales = 'free',
	strip = strip_themed(
		background_x = elem_list_rect(
			fill = rep(sp.cols, each = 2)
		) 
	)
)+theme(
	axis.text.y = element_blank(),
	axis.ticks.y = element_blank(),
	legend.key.width = unit(0.5, 'cm'),
	legend.text = element_text(size = 10),
	strip.text = element_markdown(
		size = 12,
		face = 'bold.italic',
		margin = margin(t = 10),
		lineheight = 1.1
	) 
)+scale_fill_manual(values = randomcoloR::distinctColorPalette(
	k = nrow(NU.family$summary)+1
	)
)+guides(fill = guide_legend(ncol = 1))
bar.family

ggsave(
	filename = "F2_A.jpeg",
	plot = bar.order,
	device = "jpeg",
	dpi = 300,
	width = 30,
	height = 20,
	units = 'cm'
)

ggsave(
	filename = "F2_B.jpeg",
	plot = bar.family,
	device = "jpeg",
	dpi = 300,
	width = 30,
	height = 20,
	units = 'cm'
)