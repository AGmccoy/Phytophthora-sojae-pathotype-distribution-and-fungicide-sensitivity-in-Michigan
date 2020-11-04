library(hagis)
library(readxl)
library(pander)
library(ggplot2)
library(data.table)
library(ape)
library(vegan)
library(dplyr)

citation("ape")
citation("vegan")

final_data <- read_excel("L:/Austin/!Phytophthora sojae pathotype survey/Pathotype records/Final Data set/excel Final Data set for publication.xlsx")
  
  head(final_data)
str(final_data)

  final_data$Rps <-
  
  gsub(pattern = "Rps ",
       
       replacement = "",
       
       x = final_data$Rps)
  
  final_data$Isolate <-
    gsub(pattern = "MPS17_",
         replacement = "",
         x = final_data$Isolate)

final_data$Isolate <- as.integer(final_data$Isolate)

head(final_data)

str(final_data)



hagis_args <- list(
  
  x = final_data,
  
  cutoff = 60,
  
  control = "susceptible",
  
  sample = "Isolate",
  
  gene = "Rps",
  
  perc_susc = "perc.susc"
  
)

isolates <- length(unique(final_data$Isolate))

Rps.summary <- do.call(summarize_gene, hagis_args)

pander(Rps.summary)

write.csv(Rps.summary, "rps_summary.csv", sep = ",")

autoplot(Rps.summary, type = "percentage")

autoplot(Rps.summary, type = "count")

# Figure 1 manuscript

Figure_1_data_only <- read_excel("L:/Austin/!Phytophthora sojae pathotype survey/!Pathotype paper WIP/Figures for paper/Figure 1.data_only.xlsx")

plot <- ggplot(Figure_1_data_only, aes(allele, perc.path, fill=paper))
plot <- plot + geom_bar(stat = "identity", position = 'dodge', col = '#4d4d4d')
plot <- plot + theme_gray() +
  theme(axis.text.x = element_text(size = 15, face = "bold", angle=45, hjust=1, family = "serif"),
                     axis.text.y = element_text(size = 20, face = "bold", family = "serif"),
                     axis.title.x = element_text(size = 20, face = "bold", family = "serif"),
                     axis.title.y = element_text(size = 20, face = "bold", family = "serif"),
                     axis.line.x = element_line(colour = 'gray', size=0.5, linetype='solid'),
                     axis.line.y = element_line(colour = 'gray', size=0.5, linetype='solid'),
                     legend.text = element_text(size = 10, face = "bold", family = "serif"),
                     legend.key = element_blank(),
                     legend.title = element_text(size = 10, face="bold", family = "serif"),
                     legend.position = "right",
                     strip.text.x = element_text(size = 10, face = "bold", family = "serif"),
                     title = element_text(size = 16, family = "serif")) +
                      ggtitle("Figure 1") +
  xlab("Allele ") + ylab("Percent Isolates Pathogenic")
plot


complexities <- do.call(calculate_complexities, hagis_args)

write.table(complexities$grouped_complexities, "grouped_complexities_summary.csv", sep = ",")
write.table(complexities$indvidual_complexities, "individual_complexities_summary.csv", sep = ",")


pander(complexities$grouped_complexities)



pander(complexities$indvidual_complexities)



pander(summary(complexities))



autoplot(complexities, type = "percentage")



autoplot(complexities, type = "count")


#final_data$Isolate <- as.integer(final_data$Isolate)

diversity <- do.call(calculate_diversities, hagis_args)

diversity

pander(diversity)

test <- diversities_table(diversity)                  
write.table(diversity$table_of_pathotypes, "diversity by pathotype.csv", sep = ",")
individual_pathotypes(diversity)
getwd()
ind_pathotypes <- individual_pathotypes(diversity)

write.csv(diversity$individual_pathotypes, "individual_isolate_pathotypes.csv")

##### Beta diversity analyses #####
# reading in the data
final_data <- read_excel("Final Data set/excel Final Data set for publication.xlsx", 
                                                   sheet = "kaitany and mccoy")

final_data$Rps <-
  gsub(pattern = "Rps ",
       replacement = "",
       x = final_data$Rps)


hagis_args <- list(
  
  x = final_data,
  
  cutoff = 1,
  
  control = "susceptible",
  
  sample = "Isolate",
  
  gene = "Rps",
  
  perc_susc = "Susceptible.1"
  
)

## change to binary data matrix!
final_data.matrix <- do.call(create_binary_matrix, hagis_args)

final_data.matrix


#### PCoA ####

t.diversity.pathotype.jaccard <- vegdist(final_data.matrix, "jaccard", na.rm = TRUE)


princoor.pathotype <- 
  pcoa(t.diversity.pathotype.jaccard)

## calculate the percentage of variation that each PC accounts for...
Axis1.percent <- princoor.pathotype$values$Relative_eig[[1]]*100 # Dimension (i.e., Axis 1 (PcoA1))
Axis2.percent <- princoor.pathotype$values$Relative_eig[[2]]*100 # Dimension (i.e., Axis 2 (PcoA2))

Axis1.percent
Axis2.percent

## now make a fancy looking plot that shows the PCs and the variation:
princoor.pathotype.data <- data.frame(Sample = rownames(princoor.pathotype$vectors),
                                      X = princoor.pathotype$vectors[, 1],
                                      Y = princoor.pathotype$vectors[, 2])
princoor.pathotype.data

#### adding metadata to PCoA file ####
write.csv(princoor.pathotype.data, "kaitany_mccoy_pca_data.csv")
pathotype_metadata <- read_excel("multistate_metadata.xlsx") # this file needs to be made seperately


pca.data <- left_join(princoor.pathotype.data, pathotype_metadata, by = "Sample")
write.csv(pca.data, "kaitany_mccoy_pca_data_an_metadata.csv")

#### data visualization for my survey, edited ####
ggplot(data = pca.data, aes(x = X, y = Y)) +
  geom_point(aes(colour = local)) +
  xlab(paste("PcoA1 - ", round(Axis1.percent,2), "%", sep = "")) +
  ylab(paste("PcoA2 - ", round(Axis2.percent,2), "%", sep = "")) +
  theme_bw() +
  theme(
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.key.size = unit(1, 'lines')
  ) +
  stat_ellipse(data = pca.data, aes(x = X, y = Y, color=local),
               level = 0.95) +
    ggtitle("Figure 3. Jaccard distances PCA")

# Visually it looks like there could be two pathotype groups based on study, however, we do not know 
## if these are significant yet, so we need to do a PERMANOVA.

###################################################################
# statistical tests for beta diversity                            #
# PERMANOVA, Beta-dispersion                                      #
###################################################################                                                                 

# Beta-dispersion

t.diversity.pathotype.jaccard # Jaccard distance matrix calculated above

groups <- factor(c(rep("Michigan_17",83), rep("Michigan_01", 78)))
# this makes a list of the locations for each pathotype, unsure of how it is matched up with pathotypes using betadisper...

length(groups)
length(unique(final_data$Isolate)) ## these numbers should match!!

pathotype.disp <- betadisper(t.diversity.pathotype.jaccard, groups)
anova(pathotype.disp) # tests if centroid distances are significantly different from each other
TukeyHSD(pathotype.disp) # test significance between each group
plot(pathotype.disp, hull = FALSE, ellipse = TRUE)

pathotype.disp$distances #boxplot using distances, should look similar to tukey


# PERMANOVA

pathotype.adonis <- adonis(t.diversity.pathotype.jaccard ~ groups)
pathotype.adonis

# ANOSIM analysis

pathotype.anosim <- anosim(t.diversity.pathotype.jaccard, groups)
pathotype.anosim
