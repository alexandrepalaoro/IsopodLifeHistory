rm(list=ls())

library(lme4)
library(bbmle)
library(MASS)
library(ggplot2)
library(cowplot)

fecund<-read.csv("fecundity.csv",h=T)

fecund$size.scale<-scale(log(fecund$size.fecund))
fecund$size.log<-log(fecund$size.fecund)

m1<-glmer.nb(eggs~size.scale*species+(species|pop),data=fecund)
summary(m1)

m2<-glmer(eggs~size.scale*species+(1|pop),data=fecund,family='poisson')
summary(m2)

AICtab(m3,m4)

m3<-glm(eggs~size.scale*species,data=fecund,family=poisson)
summary(m3)
plot(m3,which=1)
anova(m3,test="Chisq")

m4<-glm.nb(eggs~size.scale*species,data=fecund)
summary(m4)

m5<-glm(eggs~size.scale*species,data=fecund,family=quasipoisson)
summary(m5)

modelPlot<-glm(eggs~size.log*species,data=fecund,family=poisson)
summary(modelPlot)

coefs<-coef(modelPlot)
 
tiff(file="numbereggs.tiff",units="mm",width=170,height=120,res=600,
  compression="lzw")
           
plot(eggs~size.log,data=fecund,bty='l',las=1,cex=1.3,
     xlab="Cephalothorax width (log)",ylab="Number of eggs",
     pch=21,bg=c("black","darkgray","white")[as.numeric(species)])


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

tiff(file="eggs-diffpop.tiff",units="mm",width=170,height=120,res=600,
     compression="lzw")
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
dev.off()

par(mfrow=c(1,1))

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


tiff(file="int+slope-eggs.tiff",units="mm",width=170,height=120,res=600,
     compression="lzw")
plot_grid(plot1,plot2,labels=c("(a)","(b)"),ncol=2,nrow=1)
dev.off()
