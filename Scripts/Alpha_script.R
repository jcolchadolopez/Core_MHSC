setwd("D:/CleanData")
FaithPD <- read.table("./faith_pd.tsv", sep = "\t", header = T)
EvennES <- read.table("./evennes.tsv", sep = "\t", header = T)
ShannON <- read.table("./shannon.tsv", sep="\t", header=T)

meta <- read.csv(
	"METACONST.csv",
	header=T,
	row.names=1
)

rownames(meta)[29] <- "sample-341"

alpha.IDS <- cbind(FaithPD = FaithPD$faith_pd,
		Evenness = EvennES$pielou_evenness,
		Shannon = ShannON$shannon_entropy)

#rownames(alpha.IDS)[c(26,38)] <- c("sample-382", "sample-283")
#alpha.IDS <- alpha.IDS[order(rownames(alpha.IDS)),]
##MODIFIED
rm.nmr <- c(4,9,12,15,22,25,32,35,38,43,48,51)
alpha.IDS <- alpha.IDS[-rm.nmr,]
alpha.IDS <- cbind(meta, alpha.IDS)
alpha.IDS$species <- as.factor(alpha.IDS$species)

library("ggplot2")
library("cowplot")
library("scales")
library("ggh4x")
skipper <- strip_themed(background_x = elem_list_rect(fill = sp.cols))

bp.EEN <- ggplot(data=alpha.IDS,
		aes(x=compartment, y=Evenness, fill=species)
	)+
  	geom_boxplot()+
  	scale_fill_manual(values = sp.cols)+
	facet_wrap2(.~species, nrow=1, strip = skipper)+
	theme(
		strip.text = element_text(size = 18, face = "bold.italic"),
		plot.title = element_text(hjust = 0.5, face="bold", size=20),
		axis.title.x = element_blank(),
		axis.title.y = element_text(face="bold", size=18),
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.text.y = element_text(size = 15),
		legend.position = "none")+
	labs(title = "Host & Compartment", y=element_text("Pielou E"), tag="A")

bp.FPD <- ggplot(data=alpha.IDS,
		aes(x=compartment, y=FaithPD, fill=species)
	)+
  	geom_boxplot()+
  	scale_fill_manual(values = sp.cols)+
	facet_wrap2(.~species, nrow=1, strip = skipper)+
	theme(
		strip.text = element_text(size = 18, face = "bold.italic"),
		plot.title = element_blank(),
		axis.title.x = element_blank(),
		axis.title.y = element_text(face="bold", size=18),
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.text.y = element_text(size = 15),
		legend.position = "none")+
	labs(title = "Host & Compartment", y=element_text("Faith PD"), tag="B")

bp.SON <- ggplot(data=alpha.IDS,
		aes(x=compartment, y=Shannon, fill=species)
	)+
  	geom_boxplot()+
#  	scale_y_continuous(trans=scales::pseudo_log_trans(base=10),
#		breaks = c(100, 500, 1000, 1500, 2000)
#	)+
  	scale_fill_manual(values = sp.cols)+
	facet_wrap2(.~species, nrow=1, strip = skipper)+
	theme(
		strip.text = element_text(size = 18, face = "bold.italic"),
		plot.title = element_blank(),
		axis.title.x = element_blank(),
		axis.title.y = element_text(face="bold", size=18),
		axis.text.y = element_text(size = 15),
		axis.text.x = element_text(size = 16.5, face = "bold"),
		legend.position = "none")+
	labs(title = "Host & Compartment", y=element_text("Shannon H"), tag="C")+
	scale_x_discrete(labels=c("BS", "MR"))

plot_grid(bp.EEN, bp.FPD, bp.SON, nrow=3)

####Statistical testing
library("agricolae")
library("fitdistrplus")
dp1 <- descdist(alpha.IDS$Evenness, discrete = F)
dp2 <- descdist(alpha.IDS$FaithPD, discrete = F)
dp3 <- descdist(alpha.IDS$Shannon, discrete = F)
dList <- list(Evenness = dp1, FaithPD = dp2, Shannon = dp3)
fitW.E <- fitdist(alpha.IDS$Evenness, "weibull")
fitN.E <- fitdist(alpha.IDS$Evenness, "norm")
fitL.E <- fitdist(alpha.IDS$Evenness, "lnorm")
fitG.E <- fitdist(alpha.IDS$Evenness, "gamma")
fitO.E <- fitdist(alpha.IDS$Evenness, "logis")
ltxt <- c("Weibull", "Normal", "Lognormal", "Gamma", "Logistic")
par(mfrow = c(2,2))
cdfcomp(list(fitW.E, fitN.E, fitL.E, fitG.E, fitO.E), legendtext = ltxt)
denscomp(list(fitW.E, fitN.E, fitL.E, fitG.E, fitO.E), legendtext = ltxt)
qqcomp(list(fitW.E, fitN.E, fitL.E, fitG.E, fitO.E), legendtext = ltxt)
ppcomp(list(fitW.E, fitN.E, fitL.E, fitG.E, fitO.E), legendtext = ltxt)
##WEIBULL{'Evenness'}

fitW.F <- fitdist(alpha.IDS$FaithPD, "weibull")
fitN.F <- fitdist(alpha.IDS$FaithPD, "norm")
fitL.F <- fitdist(alpha.IDS$FaithPD, "lnorm")
fitG.F <- fitdist(alpha.IDS$FaithPD, "gamma")
fitO.F <- fitdist(alpha.IDS$FaithPD, "logis")
ltxt <- c("Weibull", "Normal", "Lognormal", "Exponential", "Gamma", "Logistic")
par(mfrow=c(2,2))
cdfcomp(list(fitW.F, fitN.F, fitL.F, fitG.F, fitO.F), legendtext = ltxt)
denscomp(list(fitW.F, fitN.F, fitL.F, fitG.F, fitO.F), legendtext = ltxt)
qqcomp(list(fitW.F, fitN.F, fitL.F, fitG.F, fitO.F), legendtext = ltxt)
ppcomp(list(fitW.F, fitN.F, fitL.F, fitG.F, fitO.F), legendtext = ltxt)
##WEIBULL|LOGNORMAL{FaithPD}

fitW.S <- fitdist(alpha.IDS$Shannon, "weibull")
fitN.S <- fitdist(alpha.IDS$Shannon, "norm")
fitL.S <- fitdist(alpha.IDS$Shannon, "lnorm")
fitG.S <- fitdist(alpha.IDS$Shannon, "gamma")
fitO.S <- fitdist(alpha.IDS$Shannon, "logis")
par(mfrow=c(2,2))
cdfcomp(list(fitW.S, fitN.S, fitL.S, fitG.S, fitO.S), legendtext = ltxt)
denscomp(list(fitW.S, fitN.S, fitL.S, fitG.S, fitO.S), legendtext = ltxt)
qqcomp(list(fitW.S, fitN.S, fitL.S, fitG.S, fitO.S), legendtext = ltxt)
ppcomp(list(fitW.S, fitN.S, fitL.S, fitG.S, fitO.S), legendtext = ltxt)
##WEIBULL|LOGISTIC{Shannon}

EVA <- glm(data = alpha.IDS,
		Evenness~species*compartment,
		family = Gamma(link = 'log')
	)
FOA <- glm(data = alpha.IDS,
		FaithPD~species*compartment,
		family = Gamma(link = 'log')
	)
SAV<- glm(data = alpha.IDS,
		Shannon~species*compartment,
		family = Gamma(link = 'log')
	)

KW.EEN <- kruskal(y = alpha.IDS$Evenness, trt = alpha.IDS$species, p.adj = "BH", group = T)
# _ X2(3) = 3.422016; p = 0.3310181; _
KW.SON <- kruskal(y = alpha.IDS$Shannon, trt = alpha.IDS$species, p.adj = "BH", group = T)
# *** X2(3) = 3.422016; p = 0.0004514015; a = hcl, b = m
KW.FPD <- kruskal(y = alpha.IDS$FaithPD, trt = alpha.IDS$species, p.adj = "BH", group = T)
# *** X2(3) = 21.3135; p = 9.061595e-5; a = h, ab = l, b = c, c = m
KW.een <- kruskal(y = alpha.IDS$Evenness, trt = alpha.IDS$compartment, p.adj = "BH", group = T)
# ** X2(2) = 9.209579; p = 0.01000381; a = B, b = SR
KW.son <- kruskal(y = alpha.IDS$Shannon, trt = alpha.IDS$compartment, p.adj = "BH", group = T)
# _ X2(2) = 5.080987; p = 0.07882749; _
KW.fpd <- kruskal(y = alpha.IDS$FaithPD, trt = alpha.IDS$compartment, p.adj = "BH", group = T)
# _ X2(2) = 1.977504; p = 0.2256541; _

####Memory
modE <- aov(glm(
modH
modF 

###Rarefaction
rarecurve(biom,
	step = 1,
	xlab = "Sample Size",
	ylab = "ASVs",
	label = T)

raremax <- min(rowSums(biom))

