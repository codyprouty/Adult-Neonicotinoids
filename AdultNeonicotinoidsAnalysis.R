##

#Analysis for adult monarch neonicotinoid manuscript

#Contact: cprouty@ufl.edu    https://scholar.google.com/citations?user=PpeDx78AAAAJ&hl=en

##

#load packages
library(afex)
library(lme4)
library(lsmeans)
library(multcomp)
library(survival)
##

setwd("C:/Users/codyp/Desktop/WorkBackup/Research/AltizerLab/AdultLD502019/")

#Import files
Master <- read.csv("AdultNeonicMay2019.csv")
Master2 <- read.csv("AdultFlightMaster.csv")
###

#Organize Data
Matings <- subset(Master, Matings > -1)
Matings1 <- subset(Master, Sex == "M")
Master$Weight <- (Master$WeightPostCages-Master$WeightEclos)
Master$WeightEq <- ((Master$WeightPostCages-Master$WeightEclos)/ Master$WeightEclos)
levels(Master2$Neonic) <- list(Ctrl="Ctrl", CL1 = "CL1", CL5 = "CL5")
Distance <- subset(Master2, Flown == 1)
Master2$WeCh <- (Master2$Day5Pre - Master2$Day1Pre) / Master2$Day1Pre
Master2$SurvObj <- with(Master2, Surv(DeadAfter, Survived == 0))
Master2$Conc <- as.factor(Master2$Conc)
###


#Experiment 1 - Lab-reared monarchs and sub-lethal effects

#Weight post cages
modelPCW <- lmer(WeightEq ~ Amount * Neonic + Sex +(1|Lin) , data = Master)
anova(modelPCW, test = "F")
PCW <- summary(modelPCW)
write.csv(PCW$coefficients, "coefficients.csv")

###

#Final - initial weight
modelTPR <- lmer(Weight ~ Amount*Neonic+Sex + (1|Lin), data = Master)
anova(modelTPR, test = "F")
TPR <- summary(modelTPR)
write.csv(TPR$coefficients, "coefficients1.csv")
###

#Average nectar consumed
modelAA <- lmer(AverageAmount ~ Amount*Neonic + Sex + (1|Lin) + WeightEclos, data = Master)
anova(modelAA, test = "F")
AA <- summary(modelAA)
write.csv(AA$coefficients, "coefficients2.csv")

###

#Egg hatching success
EHS <- mixed(EggScore ~ Amount*Neonic + WeightEclos + (1|Lin), family = poisson(link = "log"), data = Master, method="LRT")
EHS$anova_table
EHS <- summary(EHS)
write.csv(EHS$coefficients, "coefficients3.csv")

###

#Number of days survived post experiment
SUR <- mixed(TimePostRX ~ Amount*Neonic + Sex + WeightEclos + (1|Lin), family = poisson(link = "log"), data = Master, method="LRT")
SUR$anova_table
SUR <- summary(SUR)
write.csv(SUR$coefficients, "coefficients4.csv")

###

#Number of times males mated
MM <- mixed(Matings ~ Amount * Neonic + WeightEclos + (1|Lin), family=poisson(link="log"), data = Matings1, method = "LRT")
MM$anova_table
MM <- summary(MM)
write.csv(MM$coefficients, "coefficients5.csv")


tgcWeight <- summarySE(Matings1, measurevar="Matings", groupvars=c("Amount"))

plot(Matings ~ Amount, data = Matings1)

lsm<-lsmeans (MM, list( ~ Amount))
cld(lsm)

emmeans(MM, ~1, type='response', at=list(Amount=0))
emmeans(MM, ~1, type='response', at=list(Amount=500))
###

#Number of eggs laid
EGG <- mixed(Eggs ~ Amount*Neonic + WeightEclos + (1|Lin), family = poisson(link = "log"), data = Master, method="LRT")
EGG$anova_table
EGG <- summary(EGG)
write.csv(EGG$coefficients, "coefficients6.csv")


emmeans(EGG, ~1, type='response', at=list(c(Amount=0, Neonic="CL")))
emmeans(EGG, ~1, type='response', at=list(Amount=500, Neonic="CL"))

emmeans(EGG, ~1, type='response', at=list(c(Amount=0, Neonic="IM")))
emmeans(EGG, ~1, type='response', at=list(Amount=500, Neonic="IM"))

sim_bmond <- simulateResiduals(fittedModel = EGG, n = 250)
hist(sim_bmond)
plot(sim_bmond)
###

#PCA for activity observations
Sex <- Master$Sex
Group <- Master$Group
Neonic <- Master$Neonic
Lin <- Master$Lin
Amount <- Master$Amount
Feeding <- Master$Feeding
Flying <- Master$Flying
Mating <- Master$Mating
PCA <- data.frame(Sex, Group, Neonic, Amount, Lin, Feeding, Flying, Mating)

PCA <- na.omit(PCA)
#sum(is.na(PCA))
DR1 <- prcomp(PCA[-c(1,2,3,4,5)], center = TRUE, scale = TRUE)
head(DR1$x)
#there are my PCs
PCAscores <- data.frame(Sex = PCA$Sex, Group = PCA$Group, Neonic = PCA$Neonic, Conc = PCA$Amount, Lin=PCA$Lin, PC1 = DR1$x[,1], PC2 = DR1$x[,2])

PCA1POV <- DR1$sdev^2/sum(DR1$sdev^2)

GCols <- c("#999999","#FDB366","#F67E4B","#DD3D2D","#A50026","#999999","#98CAE1", "#6EA6CD", "#4A7BB7", "#364B9A")

SShapes <- c(1,2)


#Plot the PCs
plot(x = PCAscores$PC1, y = PCAscores$PC2, type = 'n', 
     xlab = paste0('PC1  (',round(PCA1POV[1]*100, digits = 2),'%)'), 
     ylab = paste0('PC2  (',round(PCA1POV[2]*100, digits = 2),'%)'),
     main = 'Principal Component Analysis'
)

for(G in unique(PCAscores$Group)){
  for(S in unique(PCAscores$Sex)){
    
    PlotSub <- PCAscores[which(PCAscores$Group == G & PCAscores$Sex == S),]
    
    TempCol <- GCols[which(unique(PCAscores$Group) == G)]
    
    TempShape <- SShapes[which(unique(PCAscores$Sex) == S)]
    
    points(x = PlotSub$PC1, y = PlotSub$PC2, col = TempCol, pch = TempShape, lwd = 2, cex = 1.5)
    
  }
}

#significance of PC1, reported in table 2
modelPCA <- lmer(PC1 ~ Conc*Neonic + Sex + (1|Lin), data = PCAscores)
anova(modelPCA, test = "F")
PCA <- summary(modelPCA)
write.csv(PCA$coefficients, "coefficients7.csv")

###

#Experiment 2 - wild monarchs and high-dose effects

#Average amount of nectar consumed
AvgA <- lm(AVGAmount ~  Conc + Source + Sex + Day1Pre, data = Master2)
anova(AvgA, test="F")
Avga <- summary(AvgA)
write.csv(Avga$coefficients, "coefficientss.csv")

emmeans(AvgA, ~1, type='response', at=list(Conc="0"))
emmeans(AvgA, ~1, type='response', at=list(Conc="1"))
emmeans(AvgA, ~1, type='response', at=list(Conc="5"))


###

#Proportional weight change during flight
WeCh <- lm(WeCh ~ Conc + Source + Sex , data = Master2)
anova(WeCh, test="F")
Wech <- summary(WeCh)
write.csv(Wech$coefficients, "coefficients1.csv")

emmeans(WeCh, ~1, type='response', at=list(Conc="0"))
emmeans(WeCh, ~1, type='response', at=list(Conc="1"))
emmeans(WeCh, ~1, type='response', at=list(Conc="5"))


###

#Distance flown
Dist <- lm(DistanceKM ~ Conc + Source + Sex + FlightChange , data = Master2)
anova(Dist, test = "F")
DisT <- summary(Dist)
write.csv(DisT$coefficients, "coefficients2.csv")

emmeans(Dist, ~1, type='response', at=list(Conc="0"))
emmeans(Dist, ~1, type='response', at=list(Conc="1"))
emmeans(Dist, ~1, type='response', at=list(Conc="5"))


###

#Average speed in flight
Speed <- lm(KMHR ~  Conc + Source + Sex + FlightChange , data = Master2)
anova(Speed, test = "F")
SPeed <- summary(Speed)
write.csv(SPeed$coefficients, "coefficients3.csv")

emmeans(Speed, ~1, type='response', at=list(Conc="0"))
emmeans(Speed, ~1, type='response', at=list(Conc="1"))
emmeans(Speed, ~1, type='response', at=list(Conc="5"))


###

#Drop test
DT <- lm(DropTest ~  Conc + Source + Sex + Day1Pre, data = Master2)
anova(DT, test = "F")
Dt <- summary(DT)
write.csv(Dt$coefficients, "coefficients4.csv")

emmeans(DT, ~1, type='response', at=list(Conc="0"))
emmeans(DT, ~1, type='response', at=list(Conc="1"))
emmeans(DT, ~1, type='response', at=list(Conc="5"))


###

#Survival
Surv <- glm(Survived ~ Conc + Sex + Day1Pre + Source, family = binomial(link="logit"), data = Master2)
anova(Surv, test="Chisq")
SUrv <- summary(Surv)
write.csv(SUrv$coefficients, "coefficients5.csv")

emmeans(Surv, ~1, type='response', at=list(Conc="0"))
emmeans(Surv, ~1, type='response', at=list(Conc="1"))
emmeans(Surv, ~1, type='response', at=list(Conc="5"))

###

#Graph of survival
Survival <- survfit(SurvObj ~ Conc, data = Master2)
plot(Survival, col=c("black","blue","darkgreen"), lty=1, xlab="Time (days)", ylab="Proportion Survived")
legend("bottomleft", c("Control", "Clothianidin 1 ppm", "Clothianidin 5 ppm"), col=c("black", "blue", "darkgreen"), lty=1, title="Treatment")
###

