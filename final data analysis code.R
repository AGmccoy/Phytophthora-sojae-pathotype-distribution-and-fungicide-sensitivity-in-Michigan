library(hagis) # Version 3.1.0
library(readxl) # Version 1.3.1
library(pander) # Version 0.6.3
library(ggplot2) # Version 3.3.0
library(data.table) # Version 1.12.8
library(ape) # Version 5.3
library(vegan) # Version 2.5-6
library(dplyr) # Version 0.8.5
library(adegenet) # Version 2.1.2
# R version 3.5.2 used for analysis

citation("ape")
citation("vegan")
citation()

final_data <-
  read_excel("excel Final Data set for publication.xlsx")

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

Figure_1_data_only <- read_excel("Figure 1.data_only.xlsx")

plot <-
  ggplot(Figure_1_data_only, aes(allele, perc.path, fill = study))
plot <-
  plot + geom_bar(stat = "identity",
                  position = 'dodge',
                  col = '#4d4d4d')
plot <- plot + theme_gray() +
  theme(
    axis.text.x = element_text(
      size = 15,
      face = "bold",
      angle = 45,
      hjust = 1,
      family = "serif"
    ),
    axis.text.y = element_text(
      size = 20,
      face = "bold",
      family = "serif"
    ),
    axis.title.x = element_text(
      size = 20,
      face = "bold",
      family = "serif"
    ),
    axis.title.y = element_text(
      size = 20,
      face = "bold",
      family = "serif"
    ),
    axis.line.x = element_line(
      colour = 'gray',
      size = 0.5,
      linetype = 'solid'
    ),
    axis.line.y = element_line(
      colour = 'gray',
      size = 0.5,
      linetype = 'solid'
    ),
    legend.text = element_text(
      size = 10,
      face = "bold",
      family = "serif"
    ),
    legend.key = element_blank(),
    legend.title = element_text(
      size = 10,
      face = "bold",
      family = "serif"
    ),
    legend.position = "right",
    strip.text.x = element_text(
      size = 15,
      face = "bold",
      family = "serif"
    ),
    title = element_text(size = 16, family = "serif")
  ) +
  ggtitle("Figure 1") +
  xlab("Allele ") + ylab("Percent Isolates Pathogenic")
plot


complexities <- do.call(calculate_complexities, hagis_args)

write.table(complexities$grouped_complexities,
            "grouped_complexities_summary.csv",
            sep = ",")
write.table(complexities$indvidual_complexities,
            "individual_complexities_summary.csv",
            sep = ",")


pander(complexities$grouped_complexities)



pander(complexities$indvidual_complexities)



pander(summary(complexities))



autoplot(complexities, type = "percentage")



autoplot(complexities, type = "count")

#### Figure 2 code ####
# we produce the file "complexity data for figure 2.xlsx" by writing out the "complexities" to a .csv file for each study and then combined them manually into one file

compare_complexities <-
  read_excel("complexity data for figure 2.xlsx")
compare_complexities$complexity <-
  as.factor(compare_complexities$complexity)

plot2 <-
  ggplot(compare_complexities,
         aes(complexity, distribution, fill = study))
plot2 <-
  plot2 + geom_bar(stat = "identity",
                   position = 'dodge',
                   col = '#4d4d4d')
plot2 <- plot2 + theme_gray() +
  theme(
    axis.text.x = element_text(
      size = 15,
      face = "bold",
      angle = 45,
      hjust = 1,
      family = "serif"
    ),
    axis.text.y = element_text(
      size = 20,
      face = "bold",
      family = "serif"
    ),
    axis.title.x = element_text(
      size = 20,
      face = "bold",
      family = "serif"
    ),
    axis.title.y = element_text(
      size = 20,
      face = "bold",
      family = "serif"
    ),
    axis.line.x = element_line(
      colour = 'gray',
      size = 0.5,
      linetype = 'solid'
    ),
    axis.line.y = element_line(
      colour = 'gray',
      size = 0.5,
      linetype = 'solid'
    ),
    legend.text = element_text(
      size = 10,
      face = "bold",
      family = "serif"
    ),
    legend.key = element_blank(),
    legend.title = element_text(
      size = 10,
      face = "bold",
      family = "serif"
    ),
    legend.position = "right",
    strip.text.x = element_text(
      size = 15,
      face = "bold",
      family = "serif"
    ),
    title = element_text(size = 16, family = "serif")
  ) +
  xlab("Complexity") + ylab("Percentage of Isolates")
plot2

#final_data$Isolate <- as.integer(final_data$Isolate)

diversity <- do.call(calculate_diversities, hagis_args)

diversity

pander(diversity)

test <- diversities_table(diversity)
write.table(diversity$table_of_pathotypes,
            "diversity by pathotype.csv",
            sep = ",")
individual_pathotypes(diversity)
getwd()
ind_pathotypes <- individual_pathotypes(diversity)

write.csv(diversity$individual_pathotypes,
          "individual_isolate_pathotypes.csv")

##### Beta diversity analyses #####
# reading in the data
final_data <-
  read_excel("excel Final Data set for publication.xlsx",
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

t.diversity.pathotype.jaccard <-
  vegdist(final_data.matrix, "jaccard", na.rm = TRUE)


princoor.pathotype <-
  pcoa(t.diversity.pathotype.jaccard)

## calculate the percentage of variation that each PC accounts for...
Axis1.percent <-
  princoor.pathotype$values$Relative_eig[[1]] * 100 # Dimension (i.e., Axis 1 (PcoA1))
Axis2.percent <-
  princoor.pathotype$values$Relative_eig[[2]] * 100 # Dimension (i.e., Axis 2 (PcoA2))

Axis1.percent
Axis2.percent

## now make a fancy looking plot that shows the PCs and the variation:
princoor.pathotype.data <-
  data.frame(
    Sample = rownames(princoor.pathotype$vectors),
    X = princoor.pathotype$vectors[, 1],
    Y = princoor.pathotype$vectors[, 2]
  )
princoor.pathotype.data

#### adding metadata to PCoA file ####
write.csv(princoor.pathotype.data, "kaitany_mccoy_pca_data.csv")
pathotype_metadata <-
  read_excel("multistate_metadata.xlsx") # this file needs to be made seperately


pca.data <-
  left_join(princoor.pathotype.data, pathotype_metadata, by = "Sample")
write.csv(pca.data, "kaitany_mccoy_pca_data_an_metadata.csv")

#### data visualization for my survey, edited ####
plot3 <- ggplot(data = pca.data, aes(x = X, y = Y)) +
  geom_point(aes(colour = study)) +
  xlab(paste("PcoA1 - ", round(Axis1.percent, 2), "%", sep = "")) +
  ylab(paste("PcoA2 - ", round(Axis2.percent, 2), "%", sep = "")) +
  theme_bw() +
  theme(
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.key.size = unit(1, 'lines')
  ) +
  stat_ellipse(data = pca.data, aes(x = X, y = Y, color = study),
               level = 0.95)
plot3

# Visually it looks like there could be two pathotype groups based on study, however, we do not know
## if these are significant yet, so we need to do a PERMANOVA.

###################################################################
# statistical tests for beta diversity                            #
# PERMANOVA, Beta-dispersion                                      #
###################################################################

# Beta-dispersion

t.diversity.pathotype.jaccard # Jaccard distance matrix calculated above

groups <- factor(c(rep("Michigan_17", 83), rep("Michigan_01", 78)))
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

# DAPC Analysis #
final_data.matrix # data for DAPC analysis is in a binary data matrix
groups <- factor(c(rep("Michigan_17", 83), rep("Michigan_01", 78)))

set.seed(999)
# using xvalDapc to identify the number of PC to retain
pathx <-
  xvalDapc(final_data.matrix, grp = groups) # this is a cross validation step used to identify the optimal number of principal coordinates to use in the analysis.
system.time(pathx2 <-
              xvalDapc(
                final_data.matrix,
                grp = groups,
                n.pca = 2:12,
                n.rep = 1000
              )) # this is the same as the step before, however, we are centering the validation around 5 PC used (2:12) and running it for 1000 reps to confirm 5 PC is the optimal number to use
pathx2

scatter(
  pathx2$DAPC,
  cex = 2,
  legend = TRUE,
  clabel = FALSE,
  posi.leg = "bottomleft",
  scree.pca = TRUE,
  posi.pca = "topleft",
  cleg = 0.75,
  xax = 1,
  yax = 2,
  inset.solid = 1
)

dapc1 <-
  dapc(
    final_data.matrix,
    var.contrib = TRUE,
    scale = FALSE,
    n.pca = 5,
    grp = groups
  ) # ran with 1 discriminant function saved
a.score.dapc1 <- a.score(dapc1) # on average only 34% accuracy
a.score.optimum.dapc1 <- optim.a.score(dapc1)

# validating the optimum PC to save for highest a-score
dapc2 <-
  dapc(final_data.matrix,
       grp = groups,
       n.pca = 14,
       n.da = 100)
optim.a.score(dapc2) # this also identifies n=5 as the optimal PC to retain for analysis

# rerunning DAPC with 5 PC as it is the optimum for this dataset
dapc3 <-
  dapc(
    final_data.matrix,
    grp = groups,
    n.pca = 5,
    var.cotrib = TRUE
  ) # 1 discriminant function saved
summary.dapc(dapc3)
a.score.dapc3 <-
  a.score(dapc3) # on average there was a 33.7% accuracy with assigning pathotypes to a priori groups
scatter(
  dapc3,
  grp = groups,
  legend = TRUE,
  scree.pca = TRUE,
  posi.pca = "topleft",
  posi.leg = "bottomleft"
)

#identifying the important resistance genes that are different between the groups (studies)
contrib <-
  loadingplot(
    dapc3$var.contr,
    axis = 1,
    threshold = 0.10,
    lab.jitter = 1
  ) # shows the resistance genes which are different from each group. threshold of 0.10 means
contrib_genes <-
  seploc(final_data.matrix) # wont work because it is not a Gen-type format, it is a matrix


## Supplementary figure 1 data and code ##

sup.fig1.data <- read_excel("Supplementary Figure 1 data.xlsx")

ggplot(sup.fig1.data, aes(x = hrs, y = od_600)) +
  stat_summary(fun = mean,
               position = position_dodge(width = 0.95),
               geom = "point") +
  #geom_line() +
  stat_smooth(method = "loess", se = FALSE, color = "black") +
  #scale_fill_manual(values=wes_palette(n=4, name="Moonrise2")) +
  stat_summary(
    fun.data = mean_se,
    position = position_dodge(width = 0.95),
    geom = "errorbar",
    col = "#4d4d4d"
  ) +
  theme_gray() +
  facet_wrap( ~ Isolate, ncol = 6) +
  guides(fill = FALSE) +
  ## scale_fill_manual(values = col[1:2]) +
  theme(
    axis.text.x = element_text(
      size = 9,
      face = "bold",
      angle = 45,
      hjust = 1,
      family = "serif"
    ),
    axis.text.y = element_text(
      size = 20,
      face = "bold",
      family = "serif"
    ),
    axis.title.x = element_text(
      size = 20,
      face = "bold",
      family = "serif"
    ),
    axis.title.y = element_text(
      size = 20,
      face = "bold",
      family = "serif"
    ),
    axis.line.x = element_line(
      colour = 'gray',
      size = 0.5,
      linetype = 'solid'
    ),
    axis.line.y = element_line(
      colour = 'gray',
      size = 0.5,
      linetype = 'solid'
    ),
    legend.text = element_text(
      size = 10,
      face = "bold",
      family = "serif"
    ),
    legend.key = element_blank(),
    legend.title = element_text(
      size = 10,
      face = "bold",
      family = "serif"
    ),
    legend.position = "right",
    strip.text.x = element_text(
      size = 10,
      face = "bold",
      family = "serif"
    ),
    title = element_text(size = 16, family = "serif")
  ) +
  xlab("hours after inoculation") + ylab("od600")
