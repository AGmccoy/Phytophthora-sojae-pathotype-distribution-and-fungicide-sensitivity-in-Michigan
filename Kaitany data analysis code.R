library(hagis)
library(readxl)
library(pander)
library(ggplot2)

Kaitany_mi <- read_excel("Kaitany_data.xlsx")

hagis_args <- list(
  
  x = Kaitany_mi,
  
  cutoff = 1,
  
  control = "Susceptible",
  
  sample = "Isolate",
  
  gene = "Rps",
  
  perc_susc = "Susceptible.1")

Rps.summary <- do.call(summarize_gene, hagis_args)

pander(Rps.summary)

autoplot(Rps.summary, type = "percentage")

autoplot(Rps.summary, type = "count")

complexities <- do.call(calculate_complexities, hagis_args)

pander(complexities$grouped_complexities)

pander(complexities$indvidual_complexities)

pander(summary(complexities))

autoplot(complexities, type = "percentage")

autoplot(complexities, type = "count")


#final_data$Isolate <- as.integer(final_data$Isolate)

diversity <- do.call(calculate_diversities, hagis_args)

diversity

pander(diversity)

diversities_table(diversity)

