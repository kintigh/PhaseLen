setwd()   
setwd("c:/users/kintigh/Dropbox (Personal)/Github/Kintigh Figs")
#One time: Be sure to install packages
#          Tools>Install Packages then add ggplot2,plyr,gmodels,deducer
#

#Each Time:In Packages tab at lower right, check plyr, ggplot2, Psych, labdsv, ca

library("ggplot2")
library("dplyr")
library("plyr")
library("stats")

setwd()   
setwd("c:/users/kintigh/Dropbox (Personal)/Github/Phaselen/Kintigh&Cowgill Graphics")


#=======================================
# FIgure 1 Generate a set of dates
#=======================================


ndates <- 25
mode="R"   #mode=N2 normal 2 sd R=rectangular
sample_sd <- 0
target_mean <- 1100
start_date <- 1000
end_date <-1200
int_len <- end_date-start_date
date_sd <-50
sample_mean <-0
Ihat_error <- 10
Ihat <- int_len + Ihat_error + 1

while ((sample_mean != target_mean) || (abs(Ihat-int_len)>=Ihat_error)) {
  seedval <- round(runif(1,min=0,max=100000))
  seedval
# seedval <-  48586  # for Table 2a, this one works but remove to retry
  set.seed(seedval)
  
  truedate <- matrix(1:ndates)
  error <- matrix(round(rnorm(ndates,mean=0,sd=date_sd)))  # get a vector of sample  measurement errors
  
  i <- 1
  
  while (i < (ndates+1)) {
    if (mode=="N2") {  #pick a true date with from a +/-2sd truncated normal distribtion 
      truedate[i] <- round(rnorm(1,mean=target_mean,sd=date_sd)) 
      if ((truedate[i]>=(target_mean-round(int_len/2))) && (truedate[i]<=(target_mean+round(int_len/2)))) i <- i+1  
    } else {           #pick a true date from a rectangular interval
      truedate[i] <- round(runif(1,min=target_mean-round(int_len/2),max=round(target_mean+round(int_len/2))))
      # print(c(i,truedate[i]))
      i <- i+1
    }
  } 
  
  measdate <- truedate+error
  sample_mean <- round(mean(measdate))
  sample_median <- median(measdate)
  sample_sd <- round(sd(measdate))
  sample_min <- min(measdate)
  sample_max <- max(measdate)

  #calculate Ihat
  if (mode == "N2"){
    K <- 4.55
  } else {
    K <- 3.46
  }
  sO <- sample_sd
  sM <- date_sd
  if (sO<sM) Ihat<- 0 
  else  {
    sT <- (sO*sO -sM*sM)^0.5
    Ihat <- round(K*sT)
  }
  cat("Random Seed=",seedval,"  Mean=",sample_mean,"  SD=",sample_sd,"  Median=",sample_median,
      "  Max=",sample_max,"  Min=",sample_min,"  Ihat=",Ihat,"\n")
}

# Write output to file for analysis by PhaseLen
sink("Table2A.ADF")
cat (ndates,", 2 #",ndates," Dates, Interval (",mode," distribution) ", start_date,"-",end_date,"  Sigma=", date_sd,"\n")
cat("# ","Random Seed=",seedval,"  Mean=",sample_mean,"  SD=",sample_sd,"  Median=",sample_median,
    "  Max=",sample_max,"  Min=",sample_min,"  Ihat=",Ihat,"\n")
for (i in 1:ndates) cat(measdate[i],50,"\n")
sink()

dates <- cbind(1:ndates,truedate,error,measdate)
colnames(dates) <- c("Sample","True","Error","Measured")
datesdf<-data.frame(dates)
datesdf

SD<-date_sd

png("Fig1.png",width=6.5,height=4,units="in",res=256)
#label <-"30 True Dates from the Interval 950-1000 with their Measured Dates"

ggplot(datesdf, aes(Sample, Measured)) +  theme_bw() +
  geom_errorbar(aes(ymin=Measured-SD, ymax=Measured+SD, width=0.25))+
  geom_point(shape=19,size=2)+
#  geom_point(aes(y=True,x=Sample), shape=9, size=2) +
#  geom_hline(yintercept=1000, size=.8) +
#  geom_hline(yintercept=1200, size=.8) +
  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_y_continuous(limits=c(900,1300),breaks=seq(900,1300,50)) + 
  coord_flip()  #+ ggtitle(label)

dev.off()

sortdf <- datesdf[order(datesdf$Measured),]
sortdf$Sample <- c(1:25)
sortdf


png("Fig1.png",width=5,height=5,units="in",res=256)
#label <-"30 True Dates from the Interval 950-1000 with their Measured Dates"

ggplot(sortdf, aes(Sample, Measured)) +  theme_bw() +
  geom_errorbar(aes(ymin=Measured-SD, ymax=Measured+SD, width=0.25))+
  geom_point(shape=19,size=2)+
#  geom_point(aes(y=True,x=Sample), shape=9, size=2) +
  geom_hline(yintercept=1000, size=.2) +
  geom_hline(yintercept=1200, size=.2) +
  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_y_continuous(limits=c(850,1350),breaks=seq(850,1350,50)) + 
  coord_flip()  #+ ggtitle(label)

dev.off()

#===============================
# Fig 4. Plot Cumulative Proportions - Empirical and Population
#===============================
pop <- read.csv("Fig4p.csv")
emp <- read.csv("Fig4e.csv")


png("Fig4u.png",width=5,height=4,units="in",res=256)
p2 <- ggplot()
p2 + 
  geom_line(data=pop, aes(x=Year, y=Prop)) +
  geom_line(data=emp, aes(x=Year, y=Prop)) +
  theme_bw() + 
  theme(axis.text.x = element_text(size=12),  axis.text.y = element_text(size=12), 
        axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_continuous(name='Year', breaks=seq(900,1300,100)) + 
  scale_y_continuous(name='Cumulative Proportion', breaks=seq(0,1,0.2))

dev.off()


#----------------------------
# Figure 5
#----------------------------
ds <- read.csv("Fig5.csv")
ds

ln <- rbind(c(6,-12),c(6,32),c(0,32))
colnames(ln) <- c("Seq","Yr")
ln <- data.frame(ln)
summary(ds)

png("Fig7.png",width=6.5,height=4,units="in",res=256)
label <-"25 Samples, 25 Measured Dates, SD=50, 200 Year True Interval"

ggplot(ds) +  theme_bw() +
  geom_point(aes(y=dmax,x=Sequence), shape=1,size=2.6)+
  geom_point(aes(y=Dist,x=Sequence), shape=5, size=2.0) +
  geom_point(aes(y=Span,x=Sequence), shape=15,size=2)+
  geom_point(aes(y=IQR,x=Sequence), shape=0, size=2) +
  geom_point(aes(y=Ihat,x=Sequence), shape=19, size=2) +
  geom_hline(yintercept=200, size=.8) +
  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_y_continuous(limits=c(0,325),breaks=seq(0,325,25),name="Interval Estimate") + 
  scale_x_continuous(name="Sample",breaks=seq(0,25,5),minor_breaks=seq(1,24,1)) +
  geom_point(aes(y=-0,x=3), shape=1,size=2.6) + geom_text(aes(y=5,x=3,label="dmax", hjust = "left"),  size=3.5) +
  geom_point(aes(y=0,x=2), shape=5, size=2.0) + geom_text(aes(y=5,x=2,label="Distance", hjust = "left"),  size=3.5) +
  geom_point(aes(y=0,x=4), shape=15, size=2) + geom_text(aes(y=5,x=4,label="Span", hjust = "left"),  size=3.5) +
  geom_point(aes(y=0,x=1), shape=0, size=2)  + geom_text(aes(y=5,x=1,label="IQR", hjust = "left"),  size=3.5) +
  geom_point(aes(y=0,x=5), shape=19, size=2)  + geom_text(aes(y=5,x=5,label="Ihat", hjust = "left"),  size=3.5) +
# geom_line(data=ln,aes(y=Yr,x=Seq))+
  coord_flip() # + ggtitle(label)


dev.off()

#png("Fig5.png",width=6.5,height=4,units="in",res=256)
label <-"25 Samples, 25 Measured Dates, SD=50, 200 Year True Interval"
subti <-"solid dot=Ihat, circle=dmax, diamond=dist triangle=Span, square=IQR"

ggplot(ds) +  theme_bw() +
  geom_point(aes(y=dmax,x=Sequence), shape=1,size=2.6)+
#  geom_point(aes(y=Dist,x=Sequence), shape=5, size=2.2) +
  geom_point(aes(y=Span,x=Sequence), shape=2,size=2)+
#  geom_point(aes(y=IQR,x=Sequence), shape=0, size=2) +
  geom_point(aes(y=Ihat,x=Sequence), shape=19, size=2) +
  geom_hline(yintercept=200, size=.8) +
  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_y_continuous(limits=c(0,325),breaks=seq(0,325,25),name="Interval Estimate") + 
  scale_x_continuous(name="Sample",breaks=seq(0,25,5),minor_breaks=seq(1,24,1)) +
  coord_flip() + ggtitle(label,subtitle=subti)

dev.off()

#============================================
# Fig 6-8
#============================================

ds <- read.csv("Test_n5-100_i0-600mn.csv")
ds

ds$Ihat<-na_if(ds$Ihat,-1)  #may need dplyr

summary(ds[c(4:9)])

summary(ds[c(4:9)])
by(ds[c(4:9)],ds$Group,summary)


png("Fig6.png",width=6.5,height=4,units="in",res=256)
g1 <- ggplot(ds,aes(y=Maxdif, x=Group)) + labs(x="Run", y="dmax") +
  scale_y_continuous(limits=c(0,600), breaks=seq(0,600,50)) + theme_bw() 
  g1 + geom_boxplot() + coord_flip()
dev.off()

png("Fig7.png",width=6.5,height=4,units="in",res=256)
g5 <- ggplot(ds,aes(y=Ihat, x=Group, )) + labs(x="Run", y="Ihat") +
  scale_y_continuous(limits=c(0,600), breaks=seq(0,600,50)) + theme_bw() 
  g5 + geom_boxplot() + coord_flip()
  dev.off()
  
png("Fig8.png",width=6.5,height=4,units="in",res=256)
g6 <- ggplot(ds,aes(y=Mid_Date, x=Group)) + labs(x="Run", y="Mean") +
  scale_y_continuous(limits=c(975,1225), name="Mean Date", breaks=seq(975,1225,25)) + theme_bw() 
  g6 + geom_boxplot() + coord_flip()
dev.off()

*************************
* Table 4
*************************
  dmaxSD<-sd(ds$Maxdif)
  distSD<-sd(ds$Dist)
  spanSD<-sd(ds$Span)
  iqrSD <-sd(ds$IQR)
  ihatSD<-sd(ds$Ihat,na.rm=TRUE)
  meanSD<-sd(ds$Mid_Date)
  dmaxM<-mean(ds$Maxdif)
  distM<-mean(ds$Dist)
  spanM<-mean(ds$Span)
  iqrM <-mean(ds$IQR)
  ihatM<-mean(ds$Ihat,na.rm=TRUE)
  meanM<-mean(ds$Mid_Date)
summary(ds,digits=5)
sds <- rbind(c(dmaxM,distM,spanM,iqrM,ihatM,meanM),c(dmaxSD,distSD,spanSD,iqrSD,ihatSD,meanSD))
colnames(sds)<-c("dmax","Dist","Span","IQR","Ihat","Mean") 
rownames(sds)<-c("Mean","SD")
print(sds,digits=2)

summary(ds$Maxdif)
by

dmaxE<-ds$Maxdif-200
distE<-ds$Dist-200
spanE<-ds$Span-200
iqrE <-ds$IQR-200
ihatE<-ds$Ihat-200
meanE<-ds$Mid_Date-1100
run  <-ds$Group

errdf<-data.frame(dmaxE,distE,spanE,iqrE,ihatE,meanE,run)
errdf
summary(errdf[c(1:6)])
by(errdf[c(1:6)],run,summary)

ddply(errdf,~run,summarize,mean=mean(dmaxE),sd=sd(dmaxE))
ddply(errdf,~run,summarize,mean=mean(distE),sd=sd(distE))
ddply(errdf,~run,summarize,mean=mean(spanE),sd=sd(spanE))
ddply(errdf,~run,summarize,mean=mean(iqrE),sd=sd(iqrE))
ddply(errdf,~run,summarize,mean=mean(ihatE,na.rm=TRUE),sd=sd(ihatE,na.rm=TRUE))
ddply(errdf,~run,summarize,mean=mean(meanE),sd=sd(meanE))


g1 <- ggplot(errdf,aes(y=dmaxE, x=run)) + labs(x="Run", y="dmax Error") +
  scale_y_continuous(limits=c(-200,200), breaks=seq(-200,200,50)) + theme_bw() 
g1 + geom_boxplot() + coord_flip()
g2 <- ggplot(errdf,aes(y=diste, x=run)) + labs(x="Run", y="Distance Error") +
  scale_y_continuous(limits=c(-200,200), breaks=seq(-200,200,50)) + theme_bw() 
g2 + geom_boxplot() + coord_flip()
g3 <- ggplot(errdf,aes(y=spanE, x=run)) + labs(x="Run", y="Span Error") +
  scale_y_continuous(limits=c(-200,200), breaks=seq(-200,200,50)) + theme_bw() 
g3 + geom_boxplot() + coord_flip()
g4 <- ggplot(errdf,aes(y=iqrE,  x=run)) + labs(x="Run", y="IQR Error") +
  scale_y_continuous(limits=c(-200,200), breaks=seq(-200,200,50)) + theme_bw() 
g4 + geom_boxplot() + coord_flip()
g5 <- ggplot(errdf,aes(y=ihatE,  x=run)) + labs(x="Run", y="Ihat Error") + 
  scale_y_continuous(limits=c(-250,250), breaks=seq(-250,250,50)) + theme_bw() 
g5 + geom_boxplot() + coord_flip()
g6 <- ggplot(errdf,aes(y=meanE, x=run)) + labs(x="Run", y="Mean Date Error") +
  scale_y_continuous(limits=c(-200,200), breaks=seq(-200,200,50)) + theme_bw() 
g6 + geom_boxplot() + coord_flip()

#=======================================
# Fig 9 & 10
#=======================================
copan <-read.csv("Fig9-CopanR.csv")
copan <-data.frame(copan)
summary(copan)
png("Fig9.png",width=6.5,height=4,units="in",res=256)
ggplot(copan, aes(Date)) + #  ggtitle(label) +
  geom_histogram(binwidth=50,breaks=seq(400,1250,50),) + theme_bw() +
#  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_x_continuous(breaks=seq(400,1250,50),limits=c(400,1250),aes(Name="Year")) + 
  scale_y_continuous(breaks=seq(0,350,50,),name="Count") 
dev.off()

C10L2TR <-read.csv("Fig10-10L2TR.csv")
C10L2TR <-data.frame(C10L2TR)
summary(C10L2TR)
png("Fig10.png",width=6.5,height=4,units="in",res=256)
ggplot(C10L2TR, aes(Date)) + #  ggtitle(label) +
  geom_histogram(breaks=seq(550,950,25)) + theme_bw() +
  #  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_x_continuous(breaks=seq(550,950,50),limits=c(550,950),aes(Name="Year")) + 
  scale_y_continuous(breaks=seq(0,12,2),limits=c(0,12),name="Count") 
dev.off()

