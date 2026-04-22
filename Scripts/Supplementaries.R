# Supplementary Table 1

require(pseq)
sup.tab.1 <- data.frame(
	Sample = sample_names(pseq),
	Sequences = sample_sums(pseq),
	AVSs = estimate_richness(pseq, measures='Observed'),
	row.names = NULL
)
write.csv(sup.tab.1, "SupplementaryTable1.csv")

# Supplementary Table 2
library("openxlsx")
wb <- createWorkbook()

addWorksheet(wb, "BetweenSites_Order")
writeData(wb, "BetweenSites_Order", h1_order %>% arrange(name))
addWorksheet(wb, "BetweenBS_Order")
writeData(wb, "BetweenBS_Order", h3_order %>% arrange(name))
addWorksheet(wb, "BetweenRS_Order")
writeData(wb, "BetweenRS_Order", h4_order %>% arrange(name))
addWorksheet(wb, "BetweenSites_Family")
writeData(wb, "BetweenSites_Family", h1_family %>% arrange(name))
addWorksheet(wb, "BetweenBS_Family")
writeData(wb, "BetweenBS_Family", h3_family %>% arrange(name))
addWorksheet(wb, "BetweenRS_Family")
writeData(wb, "BetweenRS_Family", h4_family %>% arrange(name))
saveWorkbook(wb, "SupplementaryTable2.xlsx", overwrite = TRUE)

# Supplementary Figure 1

