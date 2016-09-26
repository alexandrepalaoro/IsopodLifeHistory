rm(list=ls())

library(survival)

sobrev<-read.table("survival.csv",h=T,sep=',')

surv.model<-survfit(Surv(time,status)~species,data=sobrev,conf.type="log-log")

summary(surv.model)
plot(surv.model,col=c("red","blue"))

tiff(file="survival.tiff",units="mm",width=170,height=120,res=600,
     compression="lzw")
plot(survfit(Surv(time,status)~species,data=sobrev),bty='l',
     xlab="Days",ylab="Probability of survival",ylim=c(0,1),
     lty=c(1,2),lwd=2,las=1)
dev.off()

cox.model<-coxph(Surv(time,status)~species,data=sobrev)
summary(cox.model)
