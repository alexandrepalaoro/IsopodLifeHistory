rm(list=ls())

##
###PACKAGES NEEDED
##

library(lme4)
library(bbmle)
library(MASS)
library(ggplot2)
library(cowplot)

##LOAD DATA

fecund<-read.csv("fecundity.csv",h=T)

######################################################################################
## Scaling and centering the continuous variables that will be used in the analyses ##
######################################################################################

fecund$size.scale<-scale(log(fecund$size.fecund))
fecund$size.log<-log(fecund$size.fecund)

###################################
##    INVESTMENT IN FECUNDITY    ##
###################################

## Since we have two populations for two species, we may add population as a random effect
## First, we will add one intercept for each population, and then add one intercept for all
## populations and not distinguish between them. Lastly, we will compare the two models with 
## AIC to check which one is more informative.

m1<-glmer(eggs~size.scale*species+(species|pop),data=fecund,family="poisson")
summary(m1)

m2<-glmer(eggs~size.scale*species+(1|pop),data=fecund,family='poisson')
summary(m2)

AICtab(m1,m2)

##    dAIC df
## m2  0   7 
## m1 10   12

## The simpler model is more informative. However, before we proceed with this analysis,
## let's check if we really need to add a random effect in the data - maybe variance is 
## not that large...

#tiff(file="eggs-diffpop.tiff",units="mm",width=170,height=120,res=600,
#     compression="lzw")
par(mfrow=c(1,2))
hist(fecund$eggs[fecund$species=="floridana"&fecund$pop=="pop1"],main="",
     xlab="Number of eggs",prob=T,ylim=c(0,0.15),bty='l')
hist(fecund$eggs[fecund$species=="floridana"&fecund$pop=="pop2"],add=T,
     col=adjustcolor(2,0.6),prob=T)
text(23.67,0.15,"(a)")

hist(fecund$eggs[fecund$species=="inflata"&fecund$pop=="pop1"],xlim=c(1,15),
     main="",xlab="Number of eggs",prob=T,ylim=c(0,0.5))
hist(fecund$eggs[fecund$species=="inflata"&fecund$pop=="pop3"],add=T,col=adjustcolor("red",0.6),
     prob=T)
text(12.13,0.49,"(b)")
#dev.off()

## Indeed it isn't. Thus, in the name of parsimony, we shall use a GLM which is a much
## simpler and more straightforward analysis

par(mfrow=c(1,1))

## I performed three types of 'poisson' regressions (i.e. poisson, quasipoisson and negative binomial)
## here and then saw which one did a better job at diminishing the residual deviance.
## The poisson GLM was the best so I proceeded with it.

m3<-glm(eggs~size.scale*species,data=fecund,family=poisson)
summary(m3)
plot(m3,which=1)
anova(m3,test="Chisq")

# Negative binomial model
m4<-glm.nb(eggs~size.scale*species,data=fecund)
summary(m4)

# Quasipoisson model
m5<-glm(eggs~size.scale*species,data=fecund,family=quasipoisson)
summary(m5)

## Parameter estimation was identical with the GLMM and GLM analyses - which further highlights
## the exclusion of the random effect:

#    GLMM Parameter estimation
#Fixed effects:
#                            Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                  2.38310    0.02819   84.53  < 2e-16 ***
#size.scale                   0.39047    0.03462   11.28  < 2e-16 ***
#speciesinflata              -0.12124    0.11680   -1.04  0.29927    
#speciespetronioi            -0.31555    0.11378   -2.77  0.00555 ** 
#size.scale:speciesinflata   -0.02783    0.08903   -0.31  0.75455    
#size.scale:speciespetronioi -0.08589    0.12232   -0.70  0.48259  

#    GLM parameter estimation
#Coefficients:
#                            Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                  2.38310    0.02819  84.533  < 2e-16 ***
#size.scale                   0.39047    0.03462  11.280  < 2e-16 ***
#speciesinflata              -0.12124    0.11680  -1.038  0.29927    
#speciespetronioi            -0.31555    0.11378  -2.773  0.00555 ** 
#size.scale:speciesinflata   -0.02783    0.08903  -0.313  0.75455    
#size.scale:speciespetronioi -0.08589    0.12232  -0.702  0.48259    



# Model used for the predict() function later on.
modelPlot<-glm(eggs~size.log*species,data=fecund,family=poisson)
summary(modelPlot)

# Since I will plot using curve() we need the coeffiecients. 
coefs<-coef(modelPlot)
 
#tiff(file="numbereggs.tiff",units="mm",width=170,height=120,res=600,
#  compression="lzw")
           
plot(eggs~size.log,data=fecund,bty='l',las=1,cex=1.3,
     xlab="Cephalothorax width (log)",ylab="Number of eggs",
     pch=21,bg=c("black","darkgray","white")[as.numeric(species)])

## The poisson model formula is similar to a linear regression, in which
##                Y = a + b*X
## The main difference is an exponential function in the right side of the equation
##                Y = exp(a + b*X)
## This occurs because the link is in a log distribution. The inverse of the log 
## is an exponential, and that is why we retrieve the exponential of the coefficients
curve(exp(coefs[1]+coefs[2]*x),add=T,
      from=min(fecund$size.log[fecund$species=="floridana"]),
      to=max(fecund$size.log[fecund$species=="floridana"]),lwd=3)
curve(exp((coefs[1]+coefs[3])+(coefs[2]+coefs[5])*x),add=T,lwd=3, col="darkgray",
      from=min(fecund$size.log[fecund$species=="inflata"]),
      to=max(fecund$size.log[fecund$species=="inflata"]))
curve(exp((coefs[1]+coefs[4])+(coefs[2]+coefs[6])*x),add=T,lwd=3,lty=2,
      from=min(fecund$size.log[fecund$species=="petronioi"]),
      to=max(fecund$size.log[fecund$species=="petronioi"]))

legend("topleft", legend=c(expression(italic("Atlantoscia floridana"),
                                      italic("Atlantoscia inflata"),
                                      italic("Atlantoscia petronioi"))),
       pch = c(21,21,21),
       pt.bg=c("black","darkgray","white"),bty='n',
       cex=1.2)
dev.off()

#To estimate the parameters more precisely, we remove the (Intercept) from the R output
m3.I<-glm(eggs~size.scale*species-1,data=fecund,family=poisson)
m3.coef<-m3.I$coefficients
m3.ci<-confint(m3.I)

intercept<-data.frame(v1=c("floridana","inflata","petronioi"),
                      v2=exp(m3.ci[2:4,1]),v3=exp(m3.coef[2:4]),v4=exp(m3.ci[2:4,2]))
names(intercept)<-c("species","ci.l","intercept","ci.up")

cil.infl<-exp(m3.ci[1,1]+m3.ci[5,1])
ciup.infl<-exp(m3.ci[1,2]+m3.ci[5,2])
sl.infl<-exp(m3.coef[1]+m3.coef[5])

cil.petro<-exp(m3.ci[1,1]+m3.ci[6,1])
ciup.petro<-exp(m3.ci[1,2]+m3.ci[6,2])
sl.petro<-exp(m3.coef[1]+m3.coef[6])

slope<-data.frame(v1=c("floridana","inflata","petronioi"),
                  v2=c(exp(m3.ci[1,1]),cil.infl,cil.petro),
                  v3=c(exp(m3.coef[1]),sl.infl,sl.petro),
                  v4=c(exp(m3.ci[1,2]),ciup.infl,ciup.petro))
names(slope)<-c("species","ci.l","slope","ci.up")

#Now that all estimates are in a data.frame, let's plot them with ggplot2

plot1<-ggplot(data = intercept, aes(x = species, y = intercept, ymin = ci.l, ymax = ci.up)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(position = position_dodge(width = 0.2), width = 0.1) +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        axis.text.x = element_text(family="Arial",color = "black",size=12),
        axis.text.y = element_text(family="Arial",color="black",size=12),
        axis.title.y = element_text(family="Arial",size=13),
        axis.title.x = element_text(family="Arial",size=13),legend.position="none") +
  ylab("Intercept") + xlab("Species") +
  scale_x_discrete(breaks=c("floridana", "inflata", "petronioi"),
                   labels=c(expression(italic("A. floridana"),
                                       italic("A. inflata"),
                                       italic("A. petronioi")))) 

plot2<- ggplot(data = slope, aes(x = species, y = slope, ymin = ci.l, ymax = ci.up)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(position = position_dodge(width = 0.2), width = 0.1) +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        axis.text.x = element_text(family="Arial",color = "black",size=12),
        axis.text.y = element_text(family="Arial",color="black",size=12),
        axis.title.y = element_text(family="Arial",size=13),
        axis.title.x = element_text(family="Arial",size=13),legend.position="none") +
  ylab("Slope") + xlab("Species") +
  scale_x_discrete(breaks=c("floridana", "inflata", "petronioi"),
                   labels=c(expression(italic("A. floridana"),
                                       italic("A. inflata"),
                                       italic("A. petronioi")))) 


#tiff(file="int+slope-eggs.tiff",units="mm",width=170,height=120,res=600,
#     compression="lzw")
plot_grid(plot1,plot2,labels=c("(a)","(b)"),ncol=2,nrow=1)
#dev.off()

#DONE :D
