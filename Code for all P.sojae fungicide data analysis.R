
#install.packages("rlang")
# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.
# Source: https://gist.github.com/stevenworthington/3178163

library(plyr)
library(dplyr)
packages <- c("tidyverse", "rlang", "drc", "lme4", "lsmeans", "plotrix", "knitr", "ggplot2", "lmtest", "lmerTest", "Rmisc", "gridExtra", "plotly", "webshot", "ggpmisc", "multcompView", "growthcurver", "ggjoy", "reshape", "ggsci", "dr4pl", "purrr", "tidyverse", "xml2", "ggpubr", "ggfortify", "growthcurve", "cowplot")
lapply(packages, library, character.only = TRUE)
#devtools::install_github("briandconnelly/growthcurve", build_vignettes = TRUE)
#library(growthcurve)


## this is the code used to analyze data for the fungicide sensitivity portion of the Michigan 2017 Phytophthora survey conducted by Austin McCoy.



library(readxl)
htfs_minus_meanblank <- input file
  

  
  
## growth curve 
ggplot(htfs_minus_meanblank,aes(x = hrs, y = od_600)) +
  stat_summary(fun=mean,position=position_dodge(width=0.95), geom="point") +
  #geom_line() +
  stat_smooth(method = "loess", se = FALSE, color = "black") +
  #scale_fill_manual(values=wes_palette(n=4, name="Moonrise2")) +
  stat_summary(fun.data = mean_se,position=position_dodge(width=0.95), geom = "errorbar", col= "#4d4d4d") +
  theme_gray() +
  facet_wrap(~Isolate, ncol = 6)+
  guides(fill = FALSE) +
  ## scale_fill_manual(values = col[1:2]) +
  theme(axis.text.x = element_text(size = 9, face = "bold", angle=45, hjust=1, family = "serif"),
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
  xlab("hours ") + ylab("od600")

#function extract blank mean
extract.BLANK.mean <- function(x) {
  blank.mean <- x %>%
    subset(Isolate == "BLANK") %>%
    select(mean.od600) %>%
    as.numeric()
  return(blank.mean)
}

#function extract blank sd
extract.BLANK.sd <- function(x) {
  sd.BLANK <- x %>%
    subset(Isolate == "BLANK") %>%
    select(sd.od600) %>%
    as.numeric()
  return(sd.BLANK)
}

Z.data <- htfs_minus_meanblank %>%
  subset(conc == 0) %>%
  group_by(Isolate, chemistry, set, hrs) %>%
  summarize(mean.od600 = mean(od_600), 
            sd.od600 = sd(od_600)) %>%
  group_by(set, chemistry, hrs) %>%
  nest() %>%
  mutate(mean.BLANK = map(data, extract.BLANK.mean)) %>%
  mutate(sd.BLANK = map(data, extract.BLANK.sd)) %>%
  unnest(mean.BLANK) %>%
  unnest(sd.BLANK) %>%
  unnest(data)

Z <- function (sdPC, sdNC, mPC, mNC) {1-((3*(sdPC + sdNC))/abs(mPC - mNC))}
Z.data$Z.factor <- Z(Z.data$sd.od600, Z.data$sd.BLANK, Z.data$mean.od600, Z.data$mean.BLANK)
Z.data$Z.factor <- ifelse(Z.data$Z.factor < 0, 0, Z.data$Z.factor)
Z.data$unique <- interaction(Z.data$Isolate, Z.data$chemistry, Z.data$hrs, Z.data$set)
Z.data$chemXhrs <- interaction(Z.data$chemistry, Z.data$hrs)
Z.data$hrs <- factor(Z.data$hrs)
ggplot(Z.data, aes(x = reorder(Isolate, Z.factor), y = Z.factor, color = as.factor(hrs), shape = chemistry)) + 
  geom_point()

write.csv(Z.data, "psojae.cleanuptry2.48hr.pyraclostrobin.factorscore.csv")

## Z' filtering, saves only good isolate data for easier data analysis ##

Z.filter <- Z.data %>%
  arrange(Isolate) %>%
  subset(chemistry == "methanol") %>% ## this needs to be changed to each solvent as you go through different chemistries (i.e. acetone, water, methanol)
  select(Isolate,hrs, Z.factor, set) %>%
  spread(key = hrs, Z.factor) %>%
  na.omit()

Z.filter$filter <- ifelse(Z.filter$`48` > 0.4, "48", "Uh oh Shaggy, this needs redone")

Z.filter$unique <- ifelse(Z.filter$filter == "24", paste(Z.filter$Isolate, "mefenoxam", "24", Z.filter$set, sep = "."), ifelse(
  Z.filter$filter == "48", paste(Z.filter$Isolate, "mefenoxam", "48", Z.filter$set, sep = "."), 
  "Fail"))

write.csv(Z.filter, "psojae.cleanuptry2.pyraclostrobin.zfactorwithtime.csv")


library(readxl)
htfs.relgrowth <- read_excel("48hr/Pyraclostrobin/psojae.cleanuptry2.pyrac.48hr.ALL.xlsx")
htfs.relgrowth <- subset(htfs.relgrowth, conc != "SHAM50") # only use for pyraclostrobin to remove SHAM only plate
htfs.relgrowth <- subset(htfs.relgrowth, Isolate != "BLANK") # remove the blanks for EC50 calculation
str(htfs.relgrowth)
htfs.relgrowth$chemistry <- as.factor(htfs.relgrowth$chemistry)
htfs.relgrowth$conc <- as.numeric(htfs.relgrowth$conc)


drm.func <- function(x) {
  try(drm(relgrowth*100 ~ conc, 
          fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), 
          data = x))
}



EC50.try <- htfs.relgrowth %>%
  group_by(Isolate, set) %>% # if you add unique here it will completely throw off your code below DONT DO IT!!
  nest() %>%
  mutate(drmod = map(data, drm.func)) %>%
  mutate(class.mod = map(drmod, class)) %>%
  unnest(class.mod) 

EC50.final <- subset(EC50.try, class.mod != "try-error")

ED.func <- function(x) {
  estimate <- try(data.frame(ED(x, type = "absolute", respLev = 50, display = FALSE)))
  EC50 <- try(estimate[[1]])
  return(EC50)
}

EC50 <- EC50.final %>%
  mutate(EC50.estimate = map(drmod, ED.func)) %>%
  unnest(EC50.estimate) %>%
  select(Isolate, set, EC50.estimate)

write.csv(EC50, "psojae.cleanuptry2.48hr.pyraclostrobin.absolute.EC50.csv")


##### merging data frames to better keep track of publishable data #####
# this code merges text files for relative and absolute EC50 estimations. It was an easy way to keep all the final data together

#### Set 1 rep 1 and rep 2, preparing data for publication ####

# rep 1

rel.eth <- psojae.set1rep1.ethaboxam.relative.EC50
abs.eth <- psojae.set1rep1.ethaboxam.absolute.EC50

rel.mef <- psojae.set1rep1.mefenoxam.relative.EC50
abs.mef <- psojae.set1rep1.mefenoxam.absolute.EC50

rel.pyrac <- psojae.set1rep1.pyraclostrobin.relative.EC50
abs.pyrac <- psojae.set1rep1.pyraclostrobin.absolute.EC50

rel.oxa <- psojae.set1rep2.oxathiapiprolin.relative.EC50
abs.oxa <- psojae.set1rep2.oxathiapiprolin.absolute.EC50


rel.abs.EC50 <- merge(rel.eth, abs.eth) 
rel.abs.EC50.mef <- merge(rel.mef, abs.mef)
rel.abs.EC50.pyrac <- merge(rel.pyrac, abs.pyrac)
rel.abs.EC50.oxa <- merge(rel.oxa, abs.oxa)

Z.factor.eth.s1r1 <- psojae_ethaboxam.rep1set1_zfactor %>%
    .[, 5:10]

Z.factor.mef.s1r1 <- psojae.set1rep1.mefenoxam.zfactorscore %>%
  .[, 5:10]

Z.factor.pyrac.s1r1 <- psojae.set1rep1.pyraclostrobin.zfactorscore %>%
  .[, 5:10]

Z.factor.oxa.s1r1 <- psojae.set1rep1.oxathiapiprolin.zfactorscore %>%
  .[, 5:10]

rel.abs.EC50.zfactor <- merge(rel.abs.EC50, Z.factor.eth.s1r1)
write.csv(rel.abs.EC50.zfactor, "psojae.set1rep1.ethaboxam.alldata.csv")

rel.abs.EC50.zfactor.mef <- merge(rel.abs.EC50.mef, Z.factor.mef.s1r1)
write.csv(rel.abs.EC50.zfactor.mef, "psojae.set1rep1.mefenoxam.alldata.csv")

rel.abs.EC50.zfactor.pyrac <- merge(rel.abs.EC50.pyrac, Z.factor.pyrac.s1r1)
write.csv(rel.abs.EC50.zfactor.pyrac, "psojae.set1rep1.pyraclostrobin.alldata.csv")

rel.abs.EC50.zfactor.oxa <- merge(rel.abs.EC50.oxa, Z.factor.oxa.s1r1)
write.csv(rel.abs.EC50.zfactor.oxa, "psojae.set1rep1.oxathiapiprolin.alldata.csv")

# set 1 rep 2

rel.eth.s1r2 <- psojae.set1rep2.ethaboxam.relative.EC50
abs.eth.s1r2 <- psojae.set1rep2.ethaboxam.absolute.EC50

rel.mef.s1r2.24 <- psojae.set1rep2.redo.24hr.mefenoxam.relative.EC50
abs.mef.s1r2.24 <- psojae.set1rep2.redo.24hr.mefenoxam.absolute.EC50

rel.mef.s1r2.48 <- psojae.set1rep2.redo.48hr.mefenoxam.relative.EC50
abs.mef.s1r2.48 <- psojae.set1rep2.redo.48hr.mefenoxam.absolute.EC50

rel.pyrac.s1r2.24 <- psojae.set1rep2.redo.24hr.pyraclostrobin.relative.EC50
abs.pyrac.s1r2.24 <- psojae.set1rep2.redo.24hr.pyraclostrobin.absolute.EC50

#rel.pyrac.s1r2.48 <- 
#abs.pyrac.s1r2.48 <-  
  
rel.oxa.s1r2.24 <- psojae.set1rep2.redo.24hr.oxathiapiprolin.relative.EC50
abs.oxa.s1r2.24 <- psojae.set1rep2.redo.24hr.oxathiapiprolin.absolute.EC50

rel.oxa.s1r2.48 <- psojae.set1rep2.redo.48hr.oxathiapiprolin.relative.EC50
abs.oxa.s1r2.48 <- psojae.set1rep2.redo.48hr.oxathiapiprolin.absolute.EC50
  

rel.abs.EC50.eth.s1r2 <- merge(rel.eth.s1r2, abs.eth.s1r2) 
rel.abs.EC50.mef.s1r2.24 <- merge(rel.mef.s1r2.24, abs.mef.s1r2.24)
rel.abs.EC50.mef.s1r2.48 <- merge(rel.mef.s1r2.48, abs.mef.s1r2.48)
rel.abs.EC50.pyrac.s1r2.24 <- merge(rel.pyrac.s1r2.24, abs.pyrac.s1r2.24)
#rel.abs.EC50.pyrac.s1r2.48 <- merge(rel.pyrac.s1r2.48, abs.pyrac.s1r2.48)
rel.abs.EC50.oxa.s1r2.24 <- merge(rel.oxa.s1r2.24, abs.oxa.s1r2.24)
rel.abs.EC50.oxa.s1r2.48 <- merge(rel.oxa.s1r2.48, abs.oxa.s1r2.48)

Z.factor.eth.s1r2 <- psojae.set1rep2.ethaboxam.zfactorscore %>%
  .[, 5:10]

Z.factor.mef.s1r2.24 <- psojae.set1rep2.redos.24hr.mefenoxam.zfactorscore %>%
  .[, 5:10]

Z.factor.mef.s1r2.48 <- psojae.set1rep2.redos.48hr.mefenoxam.zfactorscore %>%
  .[, 5:10]

Z.factor.pyrac.s1r2.24 <- psojae.set1rep2.redos.24hr.pyraclostrobin.zfactorscore %>%
  .[, 5:10]

#Z.factor.pyrac.s1r2.48 <- psojae.set1rep1.pyraclostrobin.zfactorscore %>%
#  .[, 5:10]

Z.factor.oxa.s1r2.24 <- psojae.set1rep2.redos.24hr.oxathiapiprolin.zfactorscore %>%
  .[, 5:10]

Z.factor.oxa.s1r2.48 <- psojae.set1rep2.redos.48hr.oxathiapiprolin.zfactorscore %>%
  .[, 5:10]

rel.abs.EC50.zfactor.eth <- merge(rel.abs.EC50.eth.s1r2, Z.factor.eth.s1r2)
write.csv(rel.abs.EC50.zfactor.eth, "psojae.set1rep2.ethaboxam.alldata.csv")

rel.abs.EC50.zfactor.mef.24 <- merge(rel.abs.EC50.mef.s1r2.24, Z.factor.mef.s1r2.24)
write.csv(rel.abs.EC50.zfactor.mef.24, "psojae.set1rep2.24hr.mefenoxam.alldata.csv")

rel.abs.EC50.zfactor.mef.48 <- merge(rel.abs.EC50.mef.s1r2.48, Z.factor.mef.s1r2.48)
write.csv(rel.abs.EC50.zfactor.mef.48, "psojae.set1rep2.48hr.mefenoxam.alldata.csv")

rel.abs.EC50.zfactor.pyrac.24 <- merge(rel.abs.EC50.pyrac.s1r2.24, Z.factor.pyrac.s1r2.24)
write.csv(rel.abs.EC50.zfactor.pyrac.24, "psojae.set1rep1.24hr.pyraclostrobin.alldata.csv")

#rel.abs.EC50.zfactor.pyrac.48 <- merge(rel.abs.EC50.pyrac.48, Z.factor.pyrac.s1r1.48)
#write.csv(rel.abs.EC50.zfactor.pyrac.48, "psojae.set1rep1.48hr.pyraclostrobin.alldata.csv")

rel.abs.EC50.zfactor.oxa.24 <- merge(rel.abs.EC50.oxa.s1r2.24, Z.factor.oxa.s1r2.24)
write.csv(rel.abs.EC50.zfactor.oxa.24, "psojae.set1rep1.24hr.oxathiapiprolin.alldata.csv")

rel.abs.EC50.zfactor.oxa.48 <- merge(rel.abs.EC50.oxa.s1r2.48, Z.factor.oxa.s1r2.48)
write.csv(rel.abs.EC50.zfactor.oxa.48, "psojae.set1rep2.48hr.oxathiapiprolin.alldata.csv")

########################################## Set 2 ######################################################

# Rep 1

rel.eth <- psojae.set2rep1.ethaboxam.relative.EC50
abs.eth <- psojae.set2rep1.ethaboxam.absolute.EC50

rel.mef <- psojae.set2rep1.mefenoxam.relative.EC50 
abs.mef <- psojae.set2rep1.mefenoxam.absolute.EC50

rel.pyrac <- psojae.set2rep1.pyraclostrobin.relative.EC50
abs.pyrac <- psojae.set2rep1.pyraclostrobin.absolute.EC50

rel.oxa <- psojae.set2rep1.oxathiapiprolin.relative.EC50
abs.oxa <- psojae.set2rep1.oxathiapiprolin.absolute.EC50


rel.abs.EC50.eth <- merge(rel.eth, abs.eth) 
rel.abs.EC50.mef <- merge(rel.mef, abs.mef)
rel.abs.EC50.pyrac <- merge(rel.pyrac, abs.pyrac)
rel.abs.EC50.oxa <- merge(rel.oxa, abs.oxa)

Z.factor.eth.s2r1 <- psojae.set2rep1.ehtaboxam.zfactorscore %>%
  .[, 5:10]

Z.factor.mef.s2r1 <- psojae.set2rep1.mefenoxam.zfactorscore %>%
  .[, 5:10]

Z.factor.pyrac.s2r1 <- psojae.set2rep1.Pyraclostrobin.zfactorscore %>%
  .[, 5:10]

Z.factor.oxa.s2r1 <- psojae.set2rep1.Oxathiapiprolin.zfactorscore %>%
  .[, 5:10]

rel.abs.EC50.zfactor <- merge(rel.abs.EC50.eth, Z.factor.eth.s2r1)
write.csv(rel.abs.EC50.zfactor, "psojae.set2rep1.ethaboxam.alldata.csv")

rel.abs.EC50.zfactor.mef <- merge(rel.abs.EC50.mef, Z.factor.mef.s2r1)
write.csv(rel.abs.EC50.zfactor.mef, "psojae.set2rep1.mefenoxam.alldata.csv")

rel.abs.EC50.zfactor.pyrac <- merge(rel.abs.EC50.pyrac, Z.factor.pyrac.s2r1)
write.csv(rel.abs.EC50.zfactor.pyrac, "psojae.set2rep1.pyraclostrobin.alldata.csv")

rel.abs.EC50.zfactor.oxa <- merge(rel.abs.EC50.oxa, Z.factor.oxa.s2r1)
write.csv(rel.abs.EC50.zfactor.oxa, "psojae.set2rep1.oxathiapiprolin.alldata.csv")

#### rep 2 

rel.eth <- psojae.set2rep2.redo.ethaboxam.relative.EC50
abs.eth <- psojae.set2rep2.redo.ethaboxam.absolute.EC50

rel.mef <- psojae.set2rep2.mefenoxam.relative.EC50
abs.mef <- psojae.set2rep2.mefenoxam.absolute.EC50

rel.pyrac <- psojae.set2rep2.redo.pyraclostrobin.relative.EC50
abs.pyrac <- psojae.set2rep2.redo.pyraclostrobin.absolute.EC50

rel.oxa <- psojae.set2rep2.redo.oxathiapiprolin.relative.EC50
abs.oxa <- psojae.set2rep2.redo.oxathiapiprolin.absolute.EC50


rel.abs.EC50.eth <- merge(rel.eth, abs.eth) 
rel.abs.EC50.mef <- merge(rel.mef, abs.mef)
rel.abs.EC50.pyrac <- merge(rel.pyrac, abs.pyrac)
rel.abs.EC50.oxa <- merge(rel.oxa, abs.oxa)

Z.factor.eth.s2r2 <- psojae.set2rep2.redos.ethaboxam.zfactorscore %>%
  .[, 5:10]

Z.factor.mef.s2r2 <- psojae.set2rep2.mefenoxam.zfactorscore %>%
  .[, 5:10]

Z.factor.pyrac.s2r2 <- psojae.set2rep2.redos.pyraclostrobin.zfactorscore %>%
  .[, 5:10]

Z.factor.oxa.s2r2 <- psojae.set2rep2.redos.oxathiapiprolin.zfactorscore %>%
  .[, 5:10]

rel.abs.EC50.zfactor <- merge(rel.abs.EC50.eth, Z.factor.eth.s2r2)
write.csv(rel.abs.EC50.zfactor, "psojae.set2rep2.ethaboxam.alldata.csv")

rel.abs.EC50.zfactor.mef <- merge(rel.abs.EC50.mef, Z.factor.mef.s2r2)
write.csv(rel.abs.EC50.zfactor.mef, "psojae.set2rep2.mefenoxam.alldata.csv")

rel.abs.EC50.zfactor.pyrac <- merge(rel.abs.EC50.pyrac, Z.factor.pyrac.s2r2)
write.csv(rel.abs.EC50.zfactor.pyrac, "psojae.set2rep2.pyraclostrobin.alldata.csv")

rel.abs.EC50.zfactor.oxa <- merge(rel.abs.EC50.oxa, Z.factor.oxa.s2r2)
write.csv(rel.abs.EC50.zfactor.oxa, "psojae.set2rep2.oxathiapiprolin.alldata.csv")

####################################### Set 3 ###################################

# Rep 1

rel.eth <- psojae.set3rep1.48hr.ethaboxam.relative.EC50
abs.eth <- psojae.set3rep1.48hr.ethaboxam.absolute.EC50

rel.mef <- psojae.set3rep1.48hr.mefenoxam.relative.EC50
abs.mef <- psojae.set3rep1.48hr.mefenoxam.absolute.EC50

rel.pyrac <- psojae.set3rep1.48hr.pyraclostrobin.relative.EC50
abs.pyrac <- psojae.set3rep1.48hr.pyraclostrobin.absolute.EC50

rel.oxa <- psojae.set3rep1.48hr.oxathiapiprolin.relative.EC50
abs.oxa <- psojae.set3rep1.48hr.oxathiapiprolin.absolute.EC50


rel.abs.EC50.eth <- merge(rel.eth, abs.eth) 
rel.abs.EC50.mef <- merge(rel.mef, abs.mef)
rel.abs.EC50.pyrac <- merge(rel.pyrac, abs.pyrac)
rel.abs.EC50.oxa <- merge(rel.oxa, abs.oxa)

Z.factor.eth.s2r1 <- psojae.set3rep1.48hr.ethaboxam.zfactorscore %>%
  .[, 5:10]

Z.factor.mef.s2r1 <- psojae.set3rep1.48hr.mefenoxam.zfactorscore %>%
  .[, 5:10]

Z.factor.pyrac.s2r1 <- psojae.set3rep1.48hr.pyraclostrobin.zfactorscore %>%
  .[, 5:10]

Z.factor.oxa.s2r1 <- psojae.set3rep1.48hr.oxthiapiprolin.zfactorscore %>%
  .[, 5:10]

rel.abs.EC50.zfactor <- merge(rel.abs.EC50.eth, Z.factor.eth.s2r1)
write.csv(rel.abs.EC50.zfactor, "psojae.set3rep1.ethaboxam.alldata.csv")
write.csv(rel.abs.EC50.zfactor, "psojae.set3rep1.48hr.ethaboxam.alldata.csv")

rel.abs.EC50.zfactor.mef <- merge(rel.abs.EC50.mef, Z.factor.mef.s2r1)
write.csv(rel.abs.EC50.zfactor.mef, "psojae.set3rep1.mefenoxam.alldata.csv")
write.csv(rel.abs.EC50.zfactor.mef, "psojae.set3rep1.48hr.mefenoxam.alldata.csv")

rel.abs.EC50.zfactor.pyrac <- merge(rel.abs.EC50.pyrac, Z.factor.pyrac.s2r1)
write.csv(rel.abs.EC50.zfactor.pyrac, "psojae.set3rep1.pyraclostrobin.alldata.csv")
write.csv(rel.abs.EC50.zfactor.pyrac, "psojae.set3rep1.48hr.pyraclostrobin.alldata.csv")

rel.abs.EC50.zfactor.oxa <- merge(rel.abs.EC50.oxa, Z.factor.oxa.s2r1)
write.csv(rel.abs.EC50.zfactor.oxa, "psojae.set3rep1.oxathiapiprolin.alldata.csv")
write.csv(rel.abs.EC50.zfactor.oxa, "psojae.set3rep1.48hr.oxathiapiprolin.alldata.csv")

# Rep 2

rel.eth <- psojae.set3rep2r.48hr.ethaboxam.relative.EC50
abs.eth <- psojae.set3rep2r.48hr.ethaboxam.absolute.EC50

rel.mef <- psojae.set3rep2r.48hr.mefenoxam.relative.EC50
abs.mef <- psojae.set3rep2r.48hr.mefenoxam.absolute.EC50

rel.pyrac <- psojae.set3rep2r.48hr.pyraclostrobin.relative.EC50
abs.pyrac <- psojae.set3rep2r.48hr.pyraclostrobin.absolute.EC50

rel.oxa <- psojae.set3rep2r.48hr.oxathiapiprolin.relative.EC50
abs.oxa <- psojae.set3rep2r.48hr.oxathiapiprolin.absolute.EC50


rel.abs.EC50.eth <- merge(rel.eth, abs.eth) 
rel.abs.EC50.mef <- merge(rel.mef, abs.mef)
rel.abs.EC50.pyrac <- merge(rel.pyrac, abs.pyrac)
rel.abs.EC50.oxa <- merge(rel.oxa, abs.oxa)

Z.factor.eth.s3r2 <- psojae.set3rep2r.48hr.ethaboxam.factorscore %>%
  .[, 5:10]

Z.factor.mef.s3r2 <- psojae.set3rep2r.48hr.mefenoxam.zfactorscore %>%
  .[, 5:10]

Z.factor.pyrac.s3r2 <- psojae.set3rep2r.48hr.pyraclostrobin.factorscore %>%
  .[, 5:10]

Z.factor.oxa.s3r2 <- psojae.set3rep2r.48hr.oxathiapiprolin.factorscore %>%
  .[, 5:10]

rel.abs.EC50.zfactor <- merge(rel.abs.EC50.eth, Z.factor.eth.s3r2)
write.csv(rel.abs.EC50.zfactor, "psojae.set3rep2.ethaboxam.alldata.csv")
write.csv(rel.abs.EC50.zfactor, "psojae.set3rep2.48hr.ethaboxam.alldata.csv")

rel.abs.EC50.zfactor.mef <- merge(rel.abs.EC50.mef, Z.factor.mef.s3r2)
write.csv(rel.abs.EC50.zfactor.mef, "psojae.set3rep2.mefenoxam.alldata.csv")
write.csv(rel.abs.EC50.zfactor.mef, "psojae.set3rep2.48hr.mefenoxam.alldata.csv")

rel.abs.EC50.zfactor.pyrac <- merge(rel.abs.EC50.pyrac, Z.factor.pyrac.s3r2)
write.csv(rel.abs.EC50.zfactor.pyrac, "psojae.set3rep2.pyraclostrobin.alldata.csv")
write.csv(rel.abs.EC50.zfactor.pyrac, "psojae.set3rep2.48hr.pyraclostrobin.alldata.csv")

rel.abs.EC50.zfactor.oxa <- merge(rel.abs.EC50.oxa, Z.factor.oxa.s3r2)
write.csv(rel.abs.EC50.zfactor.oxa, "psojae.set3rep2.oxathiapiprolin.alldata.csv")
write.csv(rel.abs.EC50.zfactor.oxa, "psojae.set3rep2.48hr.oxathiapiprolin.alldata.csv")

####################################### cleanup try 2 ###################################

# Rep 1

rel.eth <- psojae.cleanuptry2.48hr.ethaboxam.relative.EC50
abs.eth <- psojae.cleanuptry2.48hr.ethaboxam.absolute.EC50

rel.mef <- psojae.cleanuptry2.48hr.mefenoxam.relative.EC50
abs.mef <- psojae.cleanuptry2.48hr.mefenoxam.absolute.EC50

rel.pyrac <- psojae.cleanuptry2.48hr.pyraclostrobin.relative.EC50
abs.pyrac <- psojae.cleanuptry2.48hr.pyraclostrobin.absolute.EC50

rel.oxa <- psojae.cleanuptry2.48hr.oxathiapiprolin.relative.EC50
abs.oxa <- psojae.cleanuptry2.48hr.oxathiapiprolin.absolute.EC50


rel.abs.EC50.eth <- merge(rel.eth, abs.eth) 
rel.abs.EC50.mef <- merge(rel.mef, abs.mef)
rel.abs.EC50.pyrac <- merge(rel.pyrac, abs.pyrac)
rel.abs.EC50.oxa <- merge(rel.oxa, abs.oxa)

Z.factor.eth.cu2 <- psojae.cleanuptry2.48hr.ethaboxam.factorscore %>%
  .[, 5:10]

Z.factor.mef.cu2 <- psojae.cleanuptry2.48hr.mefenoxam.factorscore %>%
  .[, 5:10]

Z.factor.pyrac.cu2 <- psojae.cleanuptry2.48hr.pyraclostrobin.factorscore %>%
  .[, 5:10]

Z.factor.oxa.cu2 <- psojae.cleanuptry2.48hr.oxathiapiprolin.factorscore %>%
  .[, 5:10]

rel.abs.EC50.zfactor <- merge(rel.abs.EC50.eth, Z.factor.eth.cu2)
write.csv(rel.abs.EC50.zfactor, "psojae.cu2.ethaboxam.alldata.csv")
write.csv(rel.abs.EC50.zfactor, "psojae.cu2.48hr.ethaboxam.alldata.csv")

rel.abs.EC50.zfactor.mef <- merge(rel.abs.EC50.mef, Z.factor.mef.cu2)
write.csv(rel.abs.EC50.zfactor.mef, "psojae.cu2.mefenoxam.alldata.csv")
write.csv(rel.abs.EC50.zfactor.mef, "psojae.cu2.48hr.mefenoxam.alldata.csv")

rel.abs.EC50.zfactor.pyrac <- merge(rel.abs.EC50.pyrac, Z.factor.pyrac.cu2)
write.csv(rel.abs.EC50.zfactor.pyrac, "psojae.cu2.pyraclostrobin.alldata.csv")
write.csv(rel.abs.EC50.zfactor.pyrac, "psojae.cu2.48hr.pyraclostrobin.alldata.csv")

rel.abs.EC50.zfactor.oxa <- merge(rel.abs.EC50.oxa, Z.factor.oxa.cu2)
write.csv(rel.abs.EC50.zfactor.oxa, "psojae.cu2.oxathiapiprolin.alldata.csv")
write.csv(rel.abs.EC50.zfactor.oxa, "psojae.cu2.48hr.oxathiapiprolin.alldata.csv")


##### summary stats for fungicide data #####
library(readxl)
library(plotrix)
fung.data.test <- read_excel("L:/Austin/!Phytophthora sojae pathotype survey/Fungicide Sensitivity/Final Fungicide Data for publication/TEST.finaldata.set1set2set3.xlsx")

fung.data.test <- subset(fung.data.test, fung.data.test$chemistry == 'pyraclostrobin') # switch to oxathiapiprolin, mefenoxam, ethaboam or pyraclostrobin to see statistics
fung.data.test$abs.EC50.estimate <- as.numeric(fung.data.test$abs.EC50.estimate)

summary(fung.data.test$abs.EC50.estimate)
stderr(fung.data.test$abs.EC50.estimate)
std.error(fung.data.test$abs.EC50.estimate)
min(fung.data.test$abs.EC50.estimate)
max(fung.data.test$abs.EC50.estimate)
mean(fung.data.test$abs.EC50.estimate)

summary(fung.data.test$rel.EC50.estimate)
min(fung.data.test$rel.EC50.estimate)
max(fung.data.test$rel.EC50.estimate)
mean(fung.data.test$rel.EC50.estimate)

length(unique(fung.data.test$Isolate))
