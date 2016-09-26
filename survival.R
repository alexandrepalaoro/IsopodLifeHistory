rm(list=ls())

##
###PACKAGES NEEDED
##

library(survival)

##LOAD DATA

sobrev<-read.table("survival.csv",h=T,sep=',')

## First, we will use the Keplan-Meier method to estimate the survival function for both species
## In our dataset,
# "time" is how long an individual lasted until death or the end of the survey;
# "status" is if the individual died or if he was alive until the end of the survey; 1 = censored, 2 = dead
# "species" is the factor with two levels

## The function Surv() creates a vector, or column, with a survival object that can be used in survival analysis

surv.model<-survfit(Surv(time,status)~species,data=sobrev,conf.type="log-log")

summary(surv.model)

## This plots the estiamted KM curve for each species; red = A. floridana, blue = A. petronioi
plot(surv.model,col=c("red","blue"))

#tiff(file="survival.tiff",units="mm",width=170,height=120,res=600,
#     compression="lzw")
plot(survfit(Surv(time,status)~species,data=sobrev),bty='l',
     xlab="Days",ylab="Probability of survival",ylim=c(0,1),
     lty=c(1,2),lwd=2,las=1)
#dev.off()

## Performing the Cox hazard model to test if the species differ in their probability of survival
cox.model<-coxph(Surv(time,status)~species,data=sobrev)
summary(cox.model)

## No, they don't

## DONE :D
