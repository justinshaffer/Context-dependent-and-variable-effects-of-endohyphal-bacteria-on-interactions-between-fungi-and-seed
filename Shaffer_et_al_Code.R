##############################################################################################################
##############################################################################################################
## R code for performing analyses done in the study: Endohyphal bacteria mitigate negative impacts of fungi on seeds of tropical trees
## 2018
## Justin P. Shaffer
## justinparkshaffer@gmail.com
##############################################################################################################
##############################################################################################################

# Set working directory and load packages
##############################################################################################################

# Replace the directory below with yours, where the data file exists.
setwd("2016_Ch2_Seed/")

# You will need these packages;
# You will also need the package "bear" which you will need to download and install a previous version of from:https://cran.r-project.org/web/packages/bear/index.html
install.packages("tidyr")
install.packages("ggplot2")
install.packages("reshape")
install.packages("reshape2")
install.packages("car")
install.packages("multcomp")
install.packages("pscl")
install.packages("betareg")
install.packages("lmtest")
install.packages("AICcmodavg")

library(tidyr)
library(ggplot2)
library(reshape)
library(reshape2)
library(car)
library(multcomp)
library(bear)
library(pscl)
library(betareg)
library(lmtest)
library(AICcmodavg)


# Read in data, add needed extra columns, and create SE file for line and bar graphs
##############################################################################################################

# Read data
data<- read.csv("Shaffer_et_al_Data.csv")

# Create columns for factor variables of otu, ehb, replicate
data$otu.f<- as.factor(data$otu)
data$ehb.f<- as.factor(data$ehb)
data$replicate.f<- as.factor(data$replicate)

# Create columns for a numeric variable for ehb
data$ehb.n<- as.numeric(data$ehb)

# Create a column for a numeric variable for fungal ecological mode (i.e., seed vs. leaf isolation source)
data$mode.n<- as.numeric(data$mode)

# Reorder values for isolate such that "control" will be last
data$isolate <- factor(data$isolate, levels=c("PS0362a", "PS0768", "PS0772", "P0265", "P0277", "P0540", "control"))

# Create standard error (SE) tables for bar and line graphs
data.col.combined.prop.SE<- summarySE(data, measurevar="col.combined.prop", groupvars=c("isolate","ehb.f","species"))
data.germ.7wk.prop.SE <- summarySE(data, measurevar="germ.7wk.prop", groupvars=c("isolate","ehb.f","species"))
data.viable.7wk.prop.SE <- summarySE(data, measurevar="viable.7wk.prop", groupvars=c("isolate","ehb.f","species"))
##############################################END OF SECTION##################################################


# Analyses for comparing seed colonization between EHB+ vs. EHB- fungal strains, within fungal isolate, for each plant species, exluding controls
##############################################################################################################

# Subset datasets to exclude controls
subset.nocontrols<- subset(data, !(isolate == "control"))
subset.nocontrols.SE<- summarySE(subset.nocontrols, measurevar = "col.combined.prop", groupvars = c("isolate", "ehb", "mode", "otu", "species", "ehb.f", "ehb.n", "otu.f", "mode.n"))

# Subset datasets to make files for each fungal isoalte
subset.nocontrols.362<- subset(subset.nocontrols, isolate=="PS0362a")
subset.nocontrols.768<- subset(subset.nocontrols, isolate=="PS0768")
subset.nocontrols.772<- subset(subset.nocontrols, isolate=="PS0772")
subset.nocontrols.265<- subset(subset.nocontrols, isolate=="P0265")
subset.nocontrols.277<- subset(subset.nocontrols, isolate=="P0277")
subset.nocontrols.540<- subset(subset.nocontrols, isolate=="P0540")

# Subset data to make files for each fungal isolate-plant species pair
subset.nocontrols.362.At<- subset(subset.nocontrols.362, species=="Apeiba_tibourbou")
subset.nocontrols.362.Cl<- subset(subset.nocontrols.362, species=="Cecropia_longipes")
subset.nocontrols.362.Cp<- subset(subset.nocontrols.362, species=="Cecropia_peltata")
subset.nocontrols.362.Op<- subset(subset.nocontrols.362, species=="Ochroma_pyramidale")
subset.nocontrols.362.Tm<- subset(subset.nocontrols.362, species=="Trema_micrantha_brown")
subset.nocontrols.768.At<- subset(subset.nocontrols.768, species=="Apeiba_tibourbou")
subset.nocontrols.768.Cl<- subset(subset.nocontrols.768, species=="Cecropia_longipes")
subset.nocontrols.768.Cp<- subset(subset.nocontrols.768, species=="Cecropia_peltata")
subset.nocontrols.768.Op<- subset(subset.nocontrols.768, species=="Ochroma_pyramidale")
subset.nocontrols.768.Tm<- subset(subset.nocontrols.768, species=="Trema_micrantha_brown")
subset.nocontrols.772.At<- subset(subset.nocontrols.772, species=="Apeiba_tibourbou")
subset.nocontrols.772.Cl<- subset(subset.nocontrols.772, species=="Cecropia_longipes")
subset.nocontrols.772.Cp<- subset(subset.nocontrols.772, species=="Cecropia_peltata")
subset.nocontrols.772.Op<- subset(subset.nocontrols.772, species=="Ochroma_pyramidale")
subset.nocontrols.772.Tm<- subset(subset.nocontrols.772, species=="Trema_micrantha_brown")
subset.nocontrols.265.At<- subset(subset.nocontrols.265, species=="Apeiba_tibourbou")
subset.nocontrols.265.Cl<- subset(subset.nocontrols.265, species=="Cecropia_longipes")
subset.nocontrols.265.Cp<- subset(subset.nocontrols.265, species=="Cecropia_peltata")
subset.nocontrols.265.Op<- subset(subset.nocontrols.265, species=="Ochroma_pyramidale")
subset.nocontrols.265.Tm<- subset(subset.nocontrols.265, species=="Trema_micrantha_brown")
subset.nocontrols.277.At<- subset(subset.nocontrols.277, species=="Apeiba_tibourbou")
subset.nocontrols.277.Cl<- subset(subset.nocontrols.277, species=="Cecropia_longipes")
subset.nocontrols.277.Cp<- subset(subset.nocontrols.277, species=="Cecropia_peltata")
subset.nocontrols.277.Op<- subset(subset.nocontrols.277, species=="Ochroma_pyramidale")
subset.nocontrols.277.Tm<- subset(subset.nocontrols.277, species=="Trema_micrantha_brown")
subset.nocontrols.540.At<- subset(subset.nocontrols.540, species=="Apeiba_tibourbou")
subset.nocontrols.540.Cl<- subset(subset.nocontrols.540, species=="Cecropia_longipes")
subset.nocontrols.540.Cp<- subset(subset.nocontrols.540, species=="Cecropia_peltata")
subset.nocontrols.540.Op<- subset(subset.nocontrols.540, species=="Ochroma_pyramidale")
subset.nocontrols.540.Tm<- subset(subset.nocontrols.540, species=="Trema_micrantha_brown")

# GRAPHING

ggplot(subset.nocontrols.SE, aes(ehb.f, y=col.combined.prop), ylim(0,1)) +
  facet_grid(species~isolate, switch = "x") +
  geom_bar(stat="identity", width=0.7, color="gray50") + 
  geom_errorbar(aes(ymin=col.combined.prop-se, ymax=col.combined.prop+se), width=0.25, colour="black", position=position_dodge(width=0.9)) +
  theme_bw() +
  theme(text=element_text(size=16), axis.title.y=element_text(size=20), axis.title.x=element_blank()) +
  ylab("Proportion seeds colonized")+
  theme(legend.position="non")+
  scale_y_continuous(limits=c(0,1.25), breaks=seq(0,1,0.25))

# ANALYSES

# Full analysis by isolate x ehb x species (beta-regression with logit link function)
subset.nocontrols$col.combined.prop.adj<-
  (subset.nocontrols$col.combined.prop*(length(subset.nocontrols$col.combined.prop)- 1) + (1/2))/(length(subset.nocontrols$col.combined.prop))
betareg.nocontrols.col.combined.prop.adj.full<-betareg(col.combined.prop.adj~isolate*ehb.f*species, data=subset.nocontrols, link = "logit")
summary(betareg.nocontrols.col.combined.prop.adj.full)
AICc(betareg.nocontrols.col.combined.prop.adj.full)

# Reduced analysis 01 - plant species only (beta-regression with logit link function)
betareg.nocontrols.col.combined.prop.adj.species<-betareg(col.combined.prop.adj~species, data=subset.nocontrols, link = "logit")
summary(betareg.nocontrols.col.combined.prop.adj.species)
AICc(betareg.nocontrols.col.combined.prop.adj.species)

# Reduced analysis 02 - fungal identity only (beta-regression with logit link function)
betareg.nocontrols.col.combined.prop.adj.isolate<-betareg(col.combined.prop.adj~isolate, data=subset.nocontrols, link = "logit")
summary(betareg.nocontrols.col.combined.prop.adj.isolate)
AICc(betareg.nocontrols.col.combined.prop.adj.isolate)

# Reduced analysis 03 - EHB status only (beta-regression with logit link function)
betareg.nocontrols.col.combined.prop.adj.ehb<-betareg(col.combined.prop.adj~ehb.f, data=subset.nocontrols, link = "logit")
summary(betareg.nocontrols.col.combined.prop.adj.ehb)
AICc(betareg.nocontrols.col.combined.prop.adj.ehb)

# Reduced analysis 04 - plant species x fungal identity (beta-regression with logit link function)
betareg.nocontrols.col.combined.prop.adj.isolate.species<-betareg(col.combined.prop.adj~isolate*species, data=subset.nocontrols, link = "logit")
summary(betareg.nocontrols.col.combined.prop.adj.isolate.species)
AICc(betareg.nocontrols.col.combined.prop.adj.isolate.species)

# Reduced analysis 05 - fungal identity x EHB status (beta-regression with logit link function)
betareg.nocontrols.col.combined.prop.adj.isolate.ehb<-betareg(col.combined.prop.adj~isolate*ehb.f, data=subset.nocontrols, link = "logit")
summary(betareg.nocontrols.col.combined.prop.adj.isolate.ehb)
AICc(betareg.nocontrols.col.combined.prop.adj.isolate.ehb)

# Likelihood ratio test - full model vs. reduced 01
lrtest(betareg.nocontrols.col.combined.prop.adj.full, 
       betareg.nocontrols.col.combined.prop.adj.species, 
       betareg.nocontrols.col.combined.prop.adj.isolate, 
       betareg.nocontrols.col.combined.prop.adj.ehb, 
       betareg.nocontrols.col.combined.prop.adj.isolate.species,
       betareg.nocontrols.col.combined.prop.adj.isolate.ehb)

# Likelihood ratio test - full model vs. reduced 02
lrtest(betareg.nocontrols.col.combined.prop.adj.full, betareg.nocontrols.col.combined.prop.adj.isolate.ehb)

# t-tests for comparing between EHB+ and EHB- strains of each fungal isolate for each plant species
t.test(subset.nocontrols.362.At$col.combined.prop~subset.nocontrols.362.At$ehb)
t.test(subset.nocontrols.362.Cl$col.combined.prop~subset.nocontrols.362.Cl$ehb)
t.test(subset.nocontrols.362.Cp$col.combined.prop~subset.nocontrols.362.Cp$ehb)
t.test(subset.nocontrols.362.Op$col.combined.prop~subset.nocontrols.362.Op$ehb)
t.test(subset.nocontrols.362.Tm$col.combined.prop~subset.nocontrols.362.Tm$ehb)
t.test(subset.nocontrols.768.At$col.combined.prop~subset.nocontrols.768.At$ehb)
t.test(subset.nocontrols.768.Cl$col.combined.prop~subset.nocontrols.768.Cl$ehb)
t.test(subset.nocontrols.768.Cp$col.combined.prop~subset.nocontrols.768.Cp$ehb)
t.test(subset.nocontrols.768.Op$col.combined.prop~subset.nocontrols.768.Op$ehb)
t.test(subset.nocontrols.768.Tm$col.combined.prop~subset.nocontrols.768.Tm$ehb)
t.test(subset.nocontrols.772.At$col.combined.prop~subset.nocontrols.772.At$ehb)
t.test(subset.nocontrols.772.Cl$col.combined.prop~subset.nocontrols.772.Cl$ehb)
t.test(subset.nocontrols.772.Cp$col.combined.prop~subset.nocontrols.772.Cp$ehb)
t.test(subset.nocontrols.772.Op$col.combined.prop~subset.nocontrols.772.Op$ehb)
t.test(subset.nocontrols.772.Tm$col.combined.prop~subset.nocontrols.772.Tm$ehb)
t.test(subset.nocontrols.265.At$col.combined.prop~subset.nocontrols.265.At$ehb)
t.test(subset.nocontrols.265.Cl$col.combined.prop~subset.nocontrols.265.Cl$ehb)
t.test(subset.nocontrols.265.Cp$col.combined.prop~subset.nocontrols.265.Cp$ehb)
t.test(subset.nocontrols.265.Op$col.combined.prop~subset.nocontrols.265.Op$ehb)
t.test(subset.nocontrols.265.Tm$col.combined.prop~subset.nocontrols.265.Tm$ehb)
t.test(subset.nocontrols.277.At$col.combined.prop~subset.nocontrols.277.At$ehb)
t.test(subset.nocontrols.277.Cl$col.combined.prop~subset.nocontrols.277.Cl$ehb)
t.test(subset.nocontrols.277.Cp$col.combined.prop~subset.nocontrols.277.Cp$ehb)
t.test(subset.nocontrols.277.Op$col.combined.prop~subset.nocontrols.277.Op$ehb)
t.test(subset.nocontrols.277.Tm$col.combined.prop~subset.nocontrols.277.Tm$ehb)
t.test(subset.nocontrols.540.At$col.combined.prop~subset.nocontrols.540.At$ehb)
t.test(subset.nocontrols.540.Cl$col.combined.prop~subset.nocontrols.540.Cl$ehb)
t.test(subset.nocontrols.540.Cp$col.combined.prop~subset.nocontrols.540.Cp$ehb)
t.test(subset.nocontrols.540.Op$col.combined.prop~subset.nocontrols.540.Op$ehb)
t.test(subset.nocontrols.540.Tm$col.combined.prop~subset.nocontrols.540.Tm$ehb)

# Create summary tables for colonization

# Treatments
mean.2014.colcombinedprop<-with(data, tapply(col.combined.prop, isolate:ehb.f:species, mean))
sd.2014.colcombinedprop<-with(data, tapply(col.combined.prop, isolate:ehb.f:species, sd))
n.2014.colcombinedprop<-with(data, tapply(col.combined.prop, isolate:ehb.f:species, length))
table.2014.colcombinedprop<-cbind(mean.2014.colcombinedprop, sd.2014.colcombinedprop, n.2014.colcombinedprop)
write.csv(table.2014.colcombinedprop, file = "SummaryTable.2014.colonization.csv")

# Controls
controls.2014<-subset(data, isolate == "control")
mean.2014.control.colcombinedprop<-with(controls.2014, tapply(col.combined.prop, isolate:species, mean))
sd.2014.control.colcombinedprop<-with(controls.2014, tapply(col.combined.prop, isolate:species, sd))
n.2014.control.colcombinedprop<-with(controls.2014, tapply(col.combined.prop, isolate:species, length))
table.2014.control.colcombinedprop<-cbind(mean.2014.control.colcombinedprop, sd.2014.control.colcombinedprop, n.2014.control.colcombinedprop)
write.csv(table.2014.control.colcombinedprop, file = "SummaryTable.2014.colonization_controls.csv")

##############################################END OF SECTION##################################################


# Analyses for comparing seed germination between seeds infected by EHB+ vs. EHB- strains of the same fungal isolate for each plant species
##############################################################################################################

# Exclude Ochroma pyramidale and Trema micrantha brown
subset.noOp.noTm<- subset(data, !(species == "Ochroma_pyramidale" | species == "Trema_micrantha_brown"))
subset.noOp.noTm$ehb[is.na(subset.noOp.noTm$ehb)] <- 0
subset.noOp.noTm$ehb.f<- as.factor(subset.noOp.noTm$ehb)
subset.noOp.noTm.germ7wkprop.SE<- summarySE(subset.noOp.noTm, measurevar = "germ.7wk.prop", groupvars = c("isolate", "mode", "otu", "ehb", "species", "ehb.f", "mode.n"))
subset.noOp.noTm.nocontrols<- subset(subset.noOp.noTm, !(isolate == "control"))

# Subset datasets to make files for each plant species, including controls
subset.At<- subset(subset.noOp.noTm, species == "Apeiba_tibourbou")
subset.Cl<- subset(subset.noOp.noTm, species == "Cecropia_longipes")
subset.Cp<- subset(subset.noOp.noTm, species == "Cecropia_peltata")
subset.At.nocontrols<- subset(subset.At, !(isolate == "control"))
subset.Cl.nocontrols<- subset(subset.Cl, !(isolate == "control"))
subset.Cp.nocontrols<- subset(subset.Cp, !(isolate == "control"))

# Subset datasets to make files for each fungal isolate, excluding controls
subset.noOp.noTm.362.nocontrols<- subset(subset.noOp.noTm, isolate=="PS0362a")
subset.noOp.noTm.768.nocontrols<- subset(subset.noOp.noTm, isolate=="PS0768")
subset.noOp.noTm.772.nocontrols<- subset(subset.noOp.noTm, isolate=="PS0772")
subset.noOp.noTm.265.nocontrols<- subset(subset.noOp.noTm, isolate=="P0265")
subset.noOp.noTm.277.nocontrols<- subset(subset.noOp.noTm, isolate=="P0277")
subset.noOp.noTm.540.nocontrols<- subset(subset.noOp.noTm, isolate=="P0540")

# Subset datasets to make files for each fungal isolate-plant species pair, excluding controls 
subset.362.At.nocontrols<- subset(subset.noOp.noTm.362.nocontrols, species=="Apeiba_tibourbou")
subset.362.Cl.nocontrols<- subset(subset.noOp.noTm.362.nocontrols, species=="Cecropia_longipes")
subset.362.Cp.nocontrols<- subset(subset.noOp.noTm.362.nocontrols, species=="Cecropia_peltata")
subset.768.At.nocontrols<- subset(subset.noOp.noTm.768.nocontrols, species=="Apeiba_tibourbou")
subset.768.Cl.nocontrols<- subset(subset.noOp.noTm.768.nocontrols, species=="Cecropia_longipes")
subset.768.Cp.nocontrols<- subset(subset.noOp.noTm.768.nocontrols, species=="Cecropia_peltata")
subset.772.At.nocontrols<- subset(subset.noOp.noTm.772.nocontrols, species=="Apeiba_tibourbou")
subset.772.Cl.nocontrols<- subset(subset.noOp.noTm.772.nocontrols, species=="Cecropia_longipes")
subset.772.Cp.nocontrols<- subset(subset.noOp.noTm.772.nocontrols, species=="Cecropia_peltata")
subset.265.At.nocontrols<- subset(subset.noOp.noTm.265.nocontrols, species=="Apeiba_tibourbou")
subset.265.Cl.nocontrols<- subset(subset.noOp.noTm.265.nocontrols, species=="Cecropia_longipes")
subset.265.Cp.nocontrols<- subset(subset.noOp.noTm.265.nocontrols, species=="Cecropia_peltata")
subset.277.At.nocontrols<- subset(subset.noOp.noTm.277.nocontrols, species=="Apeiba_tibourbou")
subset.277.Cl.nocontrols<- subset(subset.noOp.noTm.277.nocontrols, species=="Cecropia_longipes")
subset.277.Cp.nocontrols<- subset(subset.noOp.noTm.277.nocontrols, species=="Cecropia_peltata")
subset.540.At.nocontrols<- subset(subset.noOp.noTm.540.nocontrols, species=="Apeiba_tibourbou")
subset.540.Cl.nocontrols<- subset(subset.noOp.noTm.540.nocontrols, species=="Cecropia_longipes")
subset.540.Cp.nocontrols<- subset(subset.noOp.noTm.540.nocontrols, species=="Cecropia_peltata")

# GRAPHING

ggplot(subset.noOp.noTm.germ7wkprop.SE, aes(ehb.f, y=germ.7wk.prop)) +
  facet_grid(species~isolate, switch = "x", scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", width=0.7, color="gray50") + 
  geom_errorbar(aes(ymin=germ.7wk.prop-se, ymax=germ.7wk.prop+se), width=0.25, colour="black", position=position_dodge(width=0.9)) +
  theme_bw() +
  theme(text=element_text(size=16), axis.title.y=element_text(size=20), axis.title.x=element_blank()) +
  ylab("Proportion seeds germinated")+
  theme(legend.position="non")+
  scale_y_continuous(limits=c(0,1.25), breaks=seq(0,1,0.25))

# ANALYSES

# Full analysis excluding controls (generalized linear model with logit link function and binomial distribution)
glm.noOp.noTm.nocontrols.germtotalprop.glm<-glm(cbind(germ.7wk.count, ungerm.total.count-germ.7wk.count)~isolate*ehb.f*species, family=binomial(), data=subset.noOp.noTm.nocontrols)
summary(glm.noOp.noTm.nocontrols.germtotalprop.glm)
anova.noOp.noTm.nocontrols.germtotalprop.glm<- anova(glm.noOp.noTm.nocontrols.germtotalprop.glm, test="Chisq")
anova.noOp.noTm.nocontrols.germtotalprop.glm
pR2(glm.noOp.noTm.nocontrols.germtotalprop.glm)
AICc(glm.noOp.noTm.nocontrols.germtotalprop.glm)

# Full analysis excluding controls (generalized linear model with logit link function and quasibinomial distribution to detect overdispersion in the model above)
glm.noOp.noTm.nocontrols.germtotalprop.glm.quasi<-glm(cbind(germ.7wk.count, ungerm.total.count-germ.7wk.count)~isolate*ehb.f*species, family=quasibinomial(), data=subset.noOp.noTm.nocontrols)
summary(glm.noOp.noTm.nocontrols.germtotalprop.glm.quasi)

# Create generalized linear model for Apeiba tibourbou, for Dunnett's test
glm.At.germtotalprop.glm<- glm(cbind(germ.7wk.count, ungerm.total.count-germ.7wk.count)~isolate*ehb.f, family=binomial(), data=subset.At)
summary(glm.At.germtotalprop.glm)
anova.At.germtotalprop.glm<- anova(glm.At.germtotalprop.glm, test="Chisq")
anova.At.germtotalprop.glm
pR2(glm.At.germtotalprop.glm)
# Dunnett's tests (re-order isolates such that "control" is first and combine isolate and ehb columns
subset.At$isolate <- factor(subset.At$isolate, levels=c("control", "PS0362a", "PS0768", "PS0772", "P0265", "P0277", "P0540"))
subset.At$isolate.ehb.f <- paste(subset.At$isolate, subset.At$ehb.f)
subset.At$isolate.ehb.f <- as.factor(subset.At$isolate.ehb.f)
glm.At.germtotalprop.glm<- glm(cbind(germ.7wk.count, ungerm.total.count-germ.7wk.count)~isolate.ehb.f, family=binomial(), data=subset.At)
summary(glht(glm.At.germtotalprop.glm, linfct=mcp(isolate.ehb.f="Dunnett")), test=adjusted("BH"))

# Create generalized linear model for Cecropia longipes, for Dunnett's test
glm.Cl.germtotalprop.glm<- glm(cbind(germ.7wk.count, ungerm.total.count-germ.7wk.count)~isolate*ehb.f, family=binomial(), data=subset.Cl)
summary(glm.Cl.germtotalprop.glm)
anova.Cl.germtotalprop.glm<- anova(glm.Cl.germtotalprop.glm, test="Chisq")
anova.Cl.germtotalprop.glm
pR2(glm.Cl.germtotalprop.glm)
# Dunnett's tests (re-order isolates such that "control" is first and combine isolate and ehb columns
subset.Cl$isolate <- factor(subset.Cl$isolate, levels=c("control", "PS0362a", "PS0768", "PS0772", "P0265", "P0277", "P0540"))
subset.Cl$isolate.ehb.f <- paste(subset.Cl$isolate, subset.Cl$ehb.f)
subset.Cl$isolate.ehb.f <- as.factor(subset.Cl$isolate.ehb.f)
glm.Cl.germtotalprop.glm<- glm(cbind(germ.7wk.count, ungerm.total.count-germ.7wk.count)~isolate.ehb.f, family=binomial(), data=subset.Cl)
summary(glht(glm.Cl.germtotalprop.glm, linfct=mcp(isolate.ehb.f="Dunnett")), test=adjusted("BH"))

# Create generalized linear model for Cecropia peltata, for Dunnett's test
glm.Cp.germtotalprop.glm<- glm(cbind(germ.7wk.count, ungerm.total.count-germ.7wk.count)~isolate*ehb.f, family=binomial(), data=subset.Cp)
summary(glm.Cp.germtotalprop.glm)
anova.Cp.germtotalprop.glm<- anova(glm.Cp.germtotalprop.glm, test="Chisq")
anova.Cp.germtotalprop.glm
pR2(glm.Cp.germtotalprop.glm)
# Dunnett's tests (re-order isolates such that "control" is first and combine isolate and ehb columns
subset.Cp$isolate <- factor(subset.Cp$isolate, levels=c("control", "PS0362a", "PS0768", "PS0772", "P0265", "P0277", "P0540"))
subset.Cp$isolate.ehb.f <- paste(subset.Cp$isolate, subset.Cp$ehb.f)
subset.Cp$isolate.ehb.f <- as.factor(subset.Cp$isolate.ehb.f)
glm.Cp.germtotalprop.glm<- glm(cbind(germ.7wk.count, ungerm.total.count-germ.7wk.count)~isolate.ehb.f, family=binomial(), data=subset.Cp)
summary(glht(glm.Cp.germtotalprop.glm, linfct=mcp(isolate.ehb.f="Dunnett")), test=adjusted("BH"))

# t-tests for comparing between EHB+ and EHB- strains of the same fungal isolate for each plant species
t.test(subset.362.At.nocontrols$germ.7wk.prop~subset.362.At.nocontrols$ehb)
t.test(subset.362.Cl.nocontrols$germ.7wk.prop~subset.362.Cl.nocontrols$ehb)
t.test(subset.362.Cp.nocontrols$germ.7wk.prop~subset.362.Cp.nocontrols$ehb)
t.test(subset.768.At.nocontrols$germ.7wk.prop~subset.768.At.nocontrols$ehb)
t.test(subset.768.Cl.nocontrols$germ.7wk.prop~subset.768.Cl.nocontrols$ehb)
t.test(subset.768.Cp.nocontrols$germ.7wk.prop~subset.768.Cp.nocontrols$ehb)
t.test(subset.772.At.nocontrols$germ.7wk.prop~subset.772.At.nocontrols$ehb)
t.test(subset.772.Cl.nocontrols$germ.7wk.prop~subset.772.Cl.nocontrols$ehb)
t.test(subset.772.Cp.nocontrols$germ.7wk.prop~subset.772.Cp.nocontrols$ehb)
t.test(subset.265.At.nocontrols$germ.7wk.prop~subset.265.At.nocontrols$ehb)
t.test(subset.265.Cl.nocontrols$germ.7wk.prop~subset.265.Cl.nocontrols$ehb)
t.test(subset.265.Cp.nocontrols$germ.7wk.prop~subset.265.Cp.nocontrols$ehb)
t.test(subset.277.At.nocontrols$germ.7wk.prop~subset.277.At.nocontrols$ehb)
t.test(subset.277.Cl.nocontrols$germ.7wk.prop~subset.277.Cl.nocontrols$ehb)
t.test(subset.277.Cp.nocontrols$germ.7wk.prop~subset.277.Cp.nocontrols$ehb)
t.test(subset.540.At.nocontrols$germ.7wk.prop~subset.540.At.nocontrols$ehb)
t.test(subset.540.Cl.nocontrols$germ.7wk.prop~subset.540.Cl.nocontrols$ehb)
t.test(subset.540.Cp.nocontrols$germ.7wk.prop~subset.540.Cp.nocontrols$ehb)

# Create summary tables for seed germination

# Treatments
mean.2014.germtotalprop<-with(data, tapply(germ.7wk.prop, isolate:ehb.f:species, mean))
sd.2014.germtotalprop<-with(data, tapply(germ.7wk.prop, isolate:ehb.f:species, sd))
n.2014.germtotalprop<-with(data, tapply(germ.7wk.prop, isolate:ehb.f:species, length))
table.2014.germtotalprop<-cbind(mean.2014.germtotalprop, sd.2014.germtotalprop, n.2014.germtotalprop)
write.csv(table.2014.germtotalprop, file = "SummaryTable.2014.germination.csv")

# Controls
mean.2014.controls.germtotalprop<-with(data, tapply(germ.7wk.prop, isolate:species, mean))
sd.2014.controls.germtotalprop<-with(data, tapply(germ.7wk.prop, isolate:species, sd))
n.2014.controls.germtotalprop<-with(data, tapply(germ.7wk.prop, isolate:species, length))
table.2014.controls.germtotalprop<-cbind(mean.2014.controls.germtotalprop, sd.2014.controls.germtotalprop, n.2014.controls.germtotalprop)
write.csv(table.2014.controls.germtotalprop, file = "SummaryTable.2014.germination_controls.csv")

##############################################END OF SECTION##################################################


# Analyses for comparing the viability of ungerminated seeds between seeds infected by EHB+ vs. EHB- strains within the same fungalisolate, for each plant species
##############################################################################################################

# Subset datasets to exclude Ochroma pyramidale and Trema micrantha "brown"
subset.noOp.noTm<- subset(data, !(species == "Ochroma_pyramidale" | species == "Trema_micrantha_brown"))
subset.noOp.noTm$ehb[is.na(subset.noOp.noTm$ehb)] <- 0
subset.noOp.noTm$ehb.f<- as.factor(subset.noOp.noTm$ehb)
subset.noOp.noTm.viable7wkprop.SE<- summarySE(subset.noOp.noTm, measurevar = "viable.7wk.prop", groupvars = c("isolate", "mode", "otu", "ehb", "species", "ehb.f", "mode.n"))
subset.noOp.noTm.nocontrols<- subset(subset.noOp.noTm, !(isolate == "control"))

# Subset datasets to make files for each plant species, including controls
subset.At<- subset(subset.noOp.noTm, species == "Apeiba_tibourbou")
subset.Cl<- subset(subset.noOp.noTm, species == "Cecropia_longipes")
subset.Cp<- subset(subset.noOp.noTm, species == "Cecropia_peltata")
subset.At.nocontrols<- subset(subset.At, !(isolate == "control"))
subset.Cl.nocontrols<- subset(subset.Cl, !(isolate == "control"))
subset.Cp.nocontrols<- subset(subset.Cp, !(isolate == "control"))

# Subset datasets to make files for each fungal isolate, excluding controls
subset.noOp.noTm.362.nocontrols<- subset(subset.noOp.noTm, isolate=="PS0362a")
subset.noOp.noTm.768.nocontrols<- subset(subset.noOp.noTm, isolate=="PS0768")
subset.noOp.noTm.772.nocontrols<- subset(subset.noOp.noTm, isolate=="PS0772")
subset.noOp.noTm.265.nocontrols<- subset(subset.noOp.noTm, isolate=="P0265")
subset.noOp.noTm.277.nocontrols<- subset(subset.noOp.noTm, isolate=="P0277")
subset.noOp.noTm.540.nocontrols<- subset(subset.noOp.noTm, isolate=="P0540")

# Subset datasets to make files for each fungal isolate-plant species pair, excluding controls
subset.362.At.nocontrols<- subset(subset.noOp.noTm.362.nocontrols, species=="Apeiba_tibourbou")
subset.362.Cl.nocontrols<- subset(subset.noOp.noTm.362.nocontrols, species=="Cecropia_longipes")
subset.362.Cp.nocontrols<- subset(subset.noOp.noTm.362.nocontrols, species=="Cecropia_peltata")
subset.768.At.nocontrols<- subset(subset.noOp.noTm.768.nocontrols, species=="Apeiba_tibourbou")
subset.768.Cl.nocontrols<- subset(subset.noOp.noTm.768.nocontrols, species=="Cecropia_longipes")
subset.768.Cp.nocontrols<- subset(subset.noOp.noTm.768.nocontrols, species=="Cecropia_peltata")
subset.772.At.nocontrols<- subset(subset.noOp.noTm.772.nocontrols, species=="Apeiba_tibourbou")
subset.772.Cl.nocontrols<- subset(subset.noOp.noTm.772.nocontrols, species=="Cecropia_longipes")
subset.772.Cp.nocontrols<- subset(subset.noOp.noTm.772.nocontrols, species=="Cecropia_peltata")
subset.265.At.nocontrols<- subset(subset.noOp.noTm.265.nocontrols, species=="Apeiba_tibourbou")
subset.265.Cl.nocontrols<- subset(subset.noOp.noTm.265.nocontrols, species=="Cecropia_longipes")
subset.265.Cp.nocontrols<- subset(subset.noOp.noTm.265.nocontrols, species=="Cecropia_peltata")
subset.277.At.nocontrols<- subset(subset.noOp.noTm.277.nocontrols, species=="Apeiba_tibourbou")
subset.277.Cl.nocontrols<- subset(subset.noOp.noTm.277.nocontrols, species=="Cecropia_longipes")
subset.277.Cp.nocontrols<- subset(subset.noOp.noTm.277.nocontrols, species=="Cecropia_peltata")
subset.540.At.nocontrols<- subset(subset.noOp.noTm.540.nocontrols, species=="Apeiba_tibourbou")
subset.540.Cl.nocontrols<- subset(subset.noOp.noTm.540.nocontrols, species=="Cecropia_longipes")
subset.540.Cp.nocontrols<- subset(subset.noOp.noTm.540.nocontrols, species=="Cecropia_peltata")

# GRAPHING

ggplot(subset.noOp.noTm.viable7wkprop.SE, aes(ehb.f, y=viable.7wk.prop)) +
  facet_grid(species~isolate, switch = "x", scales = "free_x", space = "free_x") +
  geom_bar(stat="identity", width=0.7, color="gray50") + 
  geom_errorbar(aes(ymin=viable.7wk.prop-se, ymax=viable.7wk.prop+se), width=0.25, colour="black", position=position_dodge(width=0.9)) +
  theme_bw() +
  theme(text=element_text(size=16), axis.title.y=element_text(size=20), axis.title.x=element_blank()) +
  ylab("Proportion non-germinated seeds viable")+
  theme(legend.position="non") +
  scale_y_continuous(limits=c(0,1.25), breaks=seq(0,1,0.25))

# ANALYSES

# Full analysis excluding controls (generalized linear model with log link function and binomial distribution)
glm.noOp.noTm.nocontrols.viabletotalprop.glm<-glm(cbind(viable.7wk.count, (viable.7wk.total-viable.7wk.count))~isolate*ehb.f*species, family=binomial(), data=subset.noOp.noTm.nocontrols)
summary(glm.noOp.noTm.nocontrols.viabletotalprop.glm)
anova.noOp.noTm.nocontrols.viabletotalprop.glm<- anova(glm.noOp.noTm.nocontrols.viabletotalprop.glm, test="Chisq")
anova.noOp.noTm.nocontrols.viabletotalprop.glm
pR2(glm.noOp.noTm.nocontrols.viabletotalprop.glm)
AICc(glm.noOp.noTm.nocontrols.viabletotalprop.glm)

# Full analysis excluding controls (generalized linear model with log link function and quasibinomial distribution for checking for overdispersion in the model above)
glm.noOp.noTm.nocontrols.viabletotalprop.glm.quasi<-glm(cbind(viable.7wk.count, (viable.7wk.total-viable.7wk.count))~isolate*ehb.f*species, family=quasibinomial(), data=subset.noOp.noTm.nocontrols)
summary(glm.noOp.noTm.nocontrols.viabletotalprop.glm.quasi)

# Create linear model for Apeiba tibourbou, for Dunnett's test (will not work in GLM below due to zero variance in controls)
lm.At.viabletotalprop.lm<- lm(viable.7wk.prop~isolate*ehb.f, data=subset.At)
summary(lm.At.viabletotalprop.lm)
anova.At.viabletotalprop.lm<- aov(lm.At.viabletotalprop.lm)
summary(anova.At.viabletotalprop.lm)
# Assess normality and equal variances which appear to be fine to proceed with Dunnett's test
plot(anova.At.viabletotalprop.lm, 1)
plot(anova.At.viabletotalprop.lm, 2)
# Dunnett's tests (re-order isolates such that "control" is first and combine isolate and ehb columns)
subset.At$isolate <- factor(subset.At$isolate, levels=c("control", "PS0362a", "PS0768", "PS0772", "P0265", "P0277", "P0540"))
subset.At$isolate.ehb.f <- paste(subset.At$isolate, subset.At$ehb.f)
subset.At$isolate.ehb.f <- as.factor(subset.At$isolate.ehb.f)
lm.At.viabletotalprop.lm<- lm(viable.7wk.prop~isolate.ehb.f, data=subset.At)
anova.At.viabletotalprop.lm<- aov(lm.At.viabletotalprop.lm)
summary(glht(anova.At.viabletotalprop.lm, linfct=mcp(isolate.ehb.f="Dunnett")), test=adjusted("BH"))

# Create generalized linear model for Apeiba tibourbou, for Dunnett's test (does not work due to zero variance in controls, use linear model above)
glm.At.viabletotalprop.glm<- glm(cbind(viable.7wk.count, (viable.7wk.total-viable.7wk.count))~isolate*ehb.f, family=binomial(), data=subset.At)
summary(glm.At.viabletotalprop.glm)
anova.At.viabletotalprop.glm<- anova(glm.At.viabletotalprop.glm, test="Chisq")
anova.At.viabletotalprop.glm
pR2(glm.At.viabletotalprop.glm)
# Dunnett's tests (re-order isolates such that "control" is first and combine isolate and ehb columns)
subset.At$isolate <- factor(subset.At$isolate, levels=c("control", "PS0362a", "PS0768", "PS0772", "P0265", "P0277", "P0540"))
subset.At$isolate.ehb.f <- paste(subset.At$isolate, subset.At$ehb.f)
subset.At$isolate.ehb.f <- as.factor(subset.At$isolate.ehb.f)
glm.At.viabletotalprop.glm<- glm(cbind(viable.7wk.count, (viable.7wk.total-viable.7wk.count))~isolate.ehb.f, family=binomial(), data=subset.At)
summary(glht(glm.At.viabletotalprop.glm, linfct=mcp(isolate.ehb.f="Dunnett")), test=adjusted("BH"))

# Create generalized linear model for Cecropia longipes, for Dunnett's test
glm.Cl.viabletotalprop.glm<- glm(cbind(viable.7wk.count, (viable.7wk.total-viable.7wk.count))~isolate*ehb.f, family=binomial(), data=subset.Cl)
summary(glm.Cl.viabletotalprop.glm)
anova.Cl.viabletotalprop.glm<- anova(glm.Cl.viabletotalprop.glm, test="Chisq")
anova.Cl.viabletotalprop.glm
pR2(glm.Cl.viabletotalprop.glm)
# Dunnett's tests (re-order isolates such that "control" is first and combine isolate and ehb columns
subset.Cl$isolate <- factor(subset.Cl$isolate, levels=c("control", "PS0362a", "PS0768", "PS0772", "P0265", "P0277", "P0540"))
subset.Cl$isolate.ehb.f <- paste(subset.Cl$isolate, subset.Cl$ehb.f)
subset.Cl$isolate.ehb.f <- as.factor(subset.Cl$isolate.ehb.f)
glm.Cl.viabletotalprop.glm<- glm(cbind(viable.7wk.count, (viable.7wk.total-viable.7wk.count))~isolate.ehb.f, family=binomial(), data=subset.Cl)
summary(glht(glm.Cl.viabletotalprop.glm, linfct=mcp(isolate.ehb.f="Dunnett")), test=adjusted("BH"))

# Create generalized linear model for Cecropia peltata, for Dunnett's test
glm.Cp.viabletotalprop.glm<- glm(cbind(viable.7wk.count, (viable.7wk.total-viable.7wk.count))~isolate*ehb.f, family=binomial(), data=subset.Cp)
summary(glm.Cp.viabletotalprop.glm)
anova.Cp.viabletotalprop.glm<- anova(glm.Cp.viabletotalprop.glm, test="Chisq")
anova.Cp.viabletotalprop.glm
pR2(glm.Cp.viabletotalprop.glm)
# Dunnett's tests (re-order isolates such that "control" is first and combine isolate and ehb columns
subset.Cp$isolate <- factor(subset.Cp$isolate, levels=c("control", "PS0362a", "PS0768", "PS0772", "P0265", "P0277", "P0540"))
subset.Cp$isolate.ehb.f <- paste(subset.Cp$isolate, subset.Cp$ehb.f)
subset.Cp$isolate.ehb.f <- as.factor(subset.Cp$isolate.ehb.f)
glm.Cp.viabletotalprop.glm<- glm(cbind(viable.7wk.count, (viable.7wk.total-viable.7wk.count))~isolate.ehb.f, family=binomial(), data=subset.Cp)
summary(glht(glm.Cp.viabletotalprop.glm, linfct=mcp(isolate.ehb.f="Dunnett")), test=adjusted("BH"))

# t-tests for comparing EHB+ vs. EHB- strains within each fungal isolate for each plant species
t.test(subset.362.At.nocontrols$viable.7wk.prop~subset.362.At.nocontrols$ehb)
t.test(subset.362.Cl.nocontrols$viable.7wk.prop~subset.362.Cl.nocontrols$ehb)
t.test(subset.362.Cp.nocontrols$viable.7wk.prop~subset.362.Cp.nocontrols$ehb)
t.test(subset.768.At.nocontrols$viable.7wk.prop~subset.768.At.nocontrols$ehb)
t.test(subset.768.Cl.nocontrols$viable.7wk.prop~subset.768.Cl.nocontrols$ehb)
t.test(subset.768.Cp.nocontrols$viable.7wk.prop~subset.768.Cp.nocontrols$ehb)
t.test(subset.772.At.nocontrols$viable.7wk.prop~subset.772.At.nocontrols$ehb)
t.test(subset.772.Cl.nocontrols$viable.7wk.prop~subset.772.Cl.nocontrols$ehb)
t.test(subset.772.Cp.nocontrols$viable.7wk.prop~subset.772.Cp.nocontrols$ehb)
t.test(subset.265.At.nocontrols$viable.7wk.prop~subset.265.At.nocontrols$ehb)
t.test(subset.265.Cl.nocontrols$viable.7wk.prop~subset.265.Cl.nocontrols$ehb)
t.test(subset.265.Cp.nocontrols$viable.7wk.prop~subset.265.Cp.nocontrols$ehb)
t.test(subset.277.At.nocontrols$viable.7wk.prop~subset.277.At.nocontrols$ehb)
t.test(subset.277.Cl.nocontrols$viable.7wk.prop~subset.277.Cl.nocontrols$ehb)
t.test(subset.277.Cp.nocontrols$viable.7wk.prop~subset.277.Cp.nocontrols$ehb)
t.test(subset.540.At.nocontrols$viable.7wk.prop~subset.540.At.nocontrols$ehb)
t.test(subset.540.Cl.nocontrols$viable.7wk.prop~subset.540.Cl.nocontrols$ehb)
t.test(subset.540.Cp.nocontrols$viable.7wk.prop~subset.540.Cp.nocontrols$ehb)

# Create summary tables for the viability of ungerminated seeds

# Treatments
mean.2014.viabletotalprop<-with(data, tapply(viable.7wk.prop, isolate:ehb.f:species, mean))
sd.2014.viabletotalprop<-with(data, tapply(viable.7wk.prop, isolate:ehb.f:species, sd))
n.2014.viabletotalprop<-with(data, tapply(viable.7wk.prop, isolate:ehb.f:species, length))
table.2014.viabletotalprop<-cbind(mean.2014.viabletotalprop, sd.2014.viabletotalprop, n.2014.viabletotalprop)
write.csv(table.2014.viabletotalprop, file = "SummaryTable.2014.viability.csv")

# Controls
mean.2014.controls.viabletotalprop<-with(data, tapply(viable.7wk.prop, isolate:species, mean))
sd.2014.controls.viabletotalprop<-with(data, tapply(viable.7wk.prop, isolate:species, sd))
n.2014.controls.viabletotalprop<-with(data, tapply(viable.7wk.prop, isolate:species, length))
table.2014.controls.viabletotalprop<-cbind(mean.2014.controls.viabletotalprop, sd.2014.controls.viabletotalprop, n.2014.controls.viabletotalprop)
write.csv(table.2014.controls.viabletotalprop, file = "SummaryTable.2014.viability_controls.csv")

##############################################END OF SECTION##################################################


###############################################END OF CODE####################################################

