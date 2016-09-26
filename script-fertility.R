rm(list=ls())

##
###PACKAGES NEEDED
##

library(MASS)
library(bbmle)
library(ggplot2)
library(cowplot)

##LOAD DATA

fert<-read.csv("fertility.csv",h=T)
names(fert)

######################################################################################
## Scaling and centering the continuous variables that will be used in the analyses ##
######################################################################################

fert$size.scale<-scale(log(fert$size.fert))
fert$size.log<-log(fert$size.fert)

## First, let's take a look if populations differ in the number of embryos/mancae

#tiff(file="embryos-diffpop.tiff",units="mm",width=170,height=120,res=600,
#     compression="lzw")

par(mfrow=c(1,2))

hist(fert$embryos[fert$species=="floridana"&fert$pop=="pop1"],main="",
     xlab="Number of embryos",prob=T,ylim=c(0,0.15))
hist(fert$embryos[fert$species=="floridana"&fert$pop=="pop2"],add=T,
     col=adjustcolor(2,0.6),prob=T)
text(29,0.15,"(a)")

hist(fert$embryos[fert$species=="inflata"&fert$pop=="pop1"],xlim=c(1,15),
     main="",xlab="Number of embryos",prob=T)
hist(fert$embryos[fert$species=="inflata"&fert$pop=="pop3"],add=T,col=adjustcolor("red",0.6),
     prob=T)
text(11.01,0.29,"(b)")

#dev.off()

# No they don't. Thus, no need to perform a GLMM. Let's stick with the regular GLM.

par(mfrow=c(1,1))

# Before we go on with the analysis, there is one species that has individuals infected by the
# bacteria Wolbachia, and ones that are not. Thus, I will remove the uninfected individuals 
# because there are fewer observations of uninfected individuals, and all other analyses were 
# performed using only the infected individuals.

fert$inter<-interaction(fert$species,fert$status)

fert<-fert[!fert$inter=='petronioi.uninfected',]

# Done. Now let's do some GLMs

m1<-glm(embryos~size.scale*species,data=fert,family='poisson')
summary(m1)

m2<-glm(embryos~size.scale*species,data=fert,family=quasipoisson)
summary(m2)
anova(m2,test="Chisq")

# Data seem overdispersed and the negative binomial is not being able to estimate one of the
# parameters (i.e. the theta), so let's stick with the quasipoisson.

#model used for the predict() function
modelPlot<-glm(embryos~size.log*species,data=fert,family=quasipoisson)

#tiff(file="embryos.tiff",units="mm",width=170,height=120,res=600,
#     compression="lzw")

plot(embryos~size.log,data=fert,bty='l',las=1,cex=1.1,
     xlab="Cephalothorax width (log)",ylab="Number of embryos",
     pch=21,bg=c("black","darkgray","white")[as.numeric(species)])

legend("topleft", legend=c(expression(italic("Atlantoscia floridana"),
                                      italic("Atlantoscia inflata"),
                                      italic("Atlantoscia petronioi"))),
       pch = c(21,21,21),
       pt.bg=c("black","darkgray","white"),bty='n',
       cex=1.2)

# Here I decided to use the predict() function instead of curve() simply because I wanted
# to diversify and lear more about both :) 
# The estimated curves are identical (not shown here), so I will stick with predict() here.


# To use the predict() function we need to create a new X-axis in the same range as the plot
# and with the same length in the factors. After, we pool them all in one data.frame we use
# predict()
cw.flo<-seq(min(fert$size.log[fert$species=="floridana"]),
                max(fert$size.log[fert$species=="floridana"]),
            length.out=100)
cw.inf<-seq(min(fert$size.log[fert$species=="inflata"]),
            max(fert$size.log[fert$species=="inflata"]),
            length.out=100)
cw.pet<-seq(min(fert$size.log[fert$species=="petronioi"]),
            max(fert$size.log[fert$species=="petronioi"]),
            length.out=100)


sp.flo<-factor(rep("floridana",length(cw.flo)))
sp.inf<-factor(rep("inflata",length(cw.inf)))
sp.pet<-factor(rep("petronioi",length(cw.pet)))

Flo<-data.frame(size.log=cw.flo,species=sp.flo)
lines(Flo$size.log,predict(modelPlot,type='response',newdata=Flo),lwd=3)

infl<-data.frame(size.log=cw.inf,species=sp.inf)
lines(infl$size.log,predict(modelPlot,type='response',newdata=infl),lwd=3,col="darkgray")

Pet<-data.frame(size.log=cw.pet,species=sp.pet)
lines(Pet$size.log,predict(modelPlot,type='response',newdata=Pet),lwd=3,lty=2)

#dev.off()

#Removing the intercept to estimate parameters with more precision

m3.I<-glm(embryos~size.scale*species-1,data=fert,family=quasipoisson)
plot(m3.I,which=1,col=fert$species)
m3.coef<-m3.I$coefficients
m3.ci<-confint(m3.I)

# Calculating intercepts and slopes and adding them to a data.frame
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

#plotting an saving
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


#tiff(file="int+slope-embryos.tiff",units="mm",width=170,height=120,res=600,
#     compression="lzw")
plot_grid(plot1,plot2,labels=c("(a)","(b)"),ncol=2,nrow=1)
#dev.off()

#### Now we move to the fertility of infected vs. uninfected individuals

rm(list=ls())
fert<-read.csv("fertility.csv",h=T)
# Removing all the other species
petro<-fert[fert$species=="petronioi",]

#Scaling and centering...
petro$size.scale<-scale(log(petro$size.fert))
petro$size.log<-log(petro$size.fert)

m5<-glm(embryos~size.scale*status,data=petro,family=poisson)
summary(m5)
anova(m5,test="Chisq")

# The regular poisson GLM did a good job, so let's stick with it.

# Model for predict()
modelPlot<-glm(embryos~size.log*status,data=petro,family=poisson)

#tiff(file="status-embryos.tiff",units="mm",width=170,height=120,res=600,
#     compression="lzw")
plot(embryos~size.log,data=petro,bty='l',las=1,cex=1.3,
     xlab="Cephalothorax width (log)",ylab="Number of embryos",
     pch=21,bg=c("black","white")[as.numeric(status)])

legend("topleft", legend=c("Infected", "Non-infected"),
       pch = c(21,21),
       pt.bg=c("black","white"),bty='n',
       cex=1.2)

cw.inf<-seq(min(petro$size.log[petro$status=="infected"]),
            max(petro$size.log[petro$status=="infected"]),
            length.out=100)
cw.uninf<-seq(min(petro$size.log[petro$status=="uninfected"]),
            max(petro$size.log[petro$status=="uninfected"]),
            length.out=100)

st.inf<-factor(rep("infected",length(cw.inf)))
st.uninf<-factor(rep("uninfected",length(cw.uninf)))

Infect<-data.frame(size.log=cw.inf,status=st.inf)
lines(Infect$size.log,predict(modelPlot,type='response',newdata=Infect),lwd=3)

Uninfect<-data.frame(size.log=cw.uninf,status=st.uninf)
lines(Uninfect$size.log,predict(modelPlot,type='response',newdata=Uninfect),lwd=3,lty=2)

#dev.off()

# Calculating parameters
m5.I<-glm(embryos~size.scale*status-1,data=petro,family=poisson)
plot(m5.I,which=1)
m5.coef<-m5.I$coefficients
m5.ci<-confint(m5.I)

intercept<-data.frame(v1=c("infected","non-infected"),
                      v2=exp(m5.ci[2:3,1]),v3=exp(m5.coef[2:3]),v4=exp(m5.ci[2:3,2]))
names(intercept)<-c("status","ci.l","intercept","ci.up")

cil.uninf<-exp(m5.ci[1,1]+m5.ci[4,1])
ciup.uninf<-exp(m5.ci[1,2]+m5.ci[4,2])
sl.uninf<-exp(m5.coef[1]+m5.coef[4])

slope<-data.frame(v1=c("infected","non-infected"),
                  v2=c(exp(m5.ci[1,1]),cil.uninf),
                  v3=c(exp(m5.coef[1]),sl.uninf),
                  v4=c(exp(m5.ci[1,2]),ciup.uninf))
names(slope)<-c("status","ci.l","slope","ci.up")

# Plotting...
plot1<-ggplot(data = intercept, aes(x = status, y = intercept, ymin = ci.l, ymax = ci.up)) +
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
  ylab("Intercept") + xlab("Status") +
  scale_x_discrete(breaks=c("infected", "non-infected"),
                   labels=c("Infected","Non-infected")) 

plot2<- ggplot(data = slope, aes(x = status, y = slope, ymin = ci.l, ymax = ci.up)) +
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
  ylab("Slope") + xlab("Status") +
  scale_x_discrete(breaks=c("infected","non-infected"),
                   labels=c("Infected","Non-infected")) 


#tiff(file="int+slope-status.tiff",units="mm",width=170,height=120,res=600,
#     compression="lzw")
plot_grid(plot1,plot2,labels=c("(a)","(b)"),ncol=2,nrow=1)
#dev.off()

# There is no difference between infected and non-infected individuals
# However, I want to check the distribution of the data to see how variable it is.

#tiff(file="diff-statusANDsize.tiff",units="mm",width=170,height=120,res=600,
#    compression="lzw")
par(mfrow=c(1,2))
hist(fert$embryos[fert$species=='petronioi'&fert$status=="uninfected"],prob=T,
     ylim=c(0,0.25),main="",xlab="Number of embryos")
hist(fert$embryos[fert$species=='petronioi'&fert$status=="infected"],add=T,prob=T,
     col=adjustcolor('red',0.6))
text(13.62,0.24,"(a)")

hist(log(fert$size.fert[fert$species=='petronioi'&fert$status=="uninfected"]),prob=T,
     main="",xlab="Cepalothorax width (log)",xlim=c(0.1,0.55))
hist(log(fert$size.fert[fert$species=='petronioi'&fert$status=="infected"]),add=T,prob=T,
     col=adjustcolor('red',0.6))
text(0.48,6.7,"(b)")
#dev.off()

# Interesting to know that although there is no significant difference, infected individuals
# are less variable in the number of embryos and size than the non-infected individuals.

# DONE :D
