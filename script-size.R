rm(list=ls())

library(ggplot2)

size<-read.csv("size.csv",h=T)

tiff(file="size-diffpop.tiff",units="mm",width=170,height=120,res=600,
     compression="lzw")

par(mfrow=c(1,2))
hist(log(size$size.all[size$species=="floridana"&size$pop=="pop1"]),main="",
     xlab="Cephalothorax width (log)",prob=T,ylim=c(0,6))
hist(log(size$size.all[size$species=="floridana"&size$pop=="pop2"]),add=T,
     col=adjustcolor(2,0.6),prob=T)
text(0.54,5.91,"(a)")

hist(log(size$size.all[size$species=="inflata"&size$pop=="pop1"]),
     main="",xlab="Cephalothorax width (log)",prob=T,ylim=c(0,8))
hist(log(size$size.all[size$species=="inflata"&size$pop=="pop3"]),add=T,col=adjustcolor("red",0.6),
     prob=T)
text(0.2,7.8,"(b)")

dev.off()

par(mfrow=c(1,1))
m1<-lm(log(size.all)~species,data=size)
plot(m1)
summary(m1)
summary.aov(m1)


tiff(file="size.tiff",units="mm",width=170,height=120,res=600,
     compression="lzw")
par(bty='l')
plot(log(size.all)~species,data=size,las=1,ylab="Cephalothorax width (log)",
     xlab="Species",axes=F)
axis(1,at=c(1,2,3),labels=c(expression(italic("Atlantoscia floridana")),
                                  expression(italic("Atlantoscia inflata")),
                                  expression(italic("Atlantoscia petronioi"))))
axis(2,las=1)
box()
dev.off()

m1.I<-lm(log(size.all)~species-1,data=size)
summary(m1.I)
m1.ci<-confint(m1.I)
m1.coef<-m1.I$coefficients

intercept<-data.frame(v1=c("floridana","inflata","petronioi"),v2=m1.ci[,1],v3=m1.coef,v4=m1.ci[,2])
names(intercept)<-c("species","ci.l","intercept","ci.up")

tiff(file="intercept-size.tiff",units="mm",width=170,height=120,res=600,
     compression="lzw")
ggplot(data = intercept, aes(x = species, y = intercept, ymin = ci.l, ymax = ci.up)) +
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
dev.off()
