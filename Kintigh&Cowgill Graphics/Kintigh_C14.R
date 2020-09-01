setwd()   
setwd("c:/users/kintigh/dropbox (ASU)/Archy Research/C14/AmAntiq/Kintigh Figs")
#One time: Be sure to install packages
#          Tools>Install Packages then add ggplot2,plyr,gmodels,deducer
#

#Each Time:In Packages tab at lower right, check plyr, ggplot2, Psych, labdsv, ca

library("ggplot2")
library("dplyr")
library("plyr")
library("stats")

setwd("c:/users/kintigh/dropbox (ASU)/Archy Research/C14/AmAntiq/Kintigh Figs")
getwd() 


#=======================================
# FIgure 1, Table 2 - Generate a set of dates
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


png("Fig1a.png",width=5,height=5,units="in",res=256)
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

#----------------------------
# Figure 7
#----------------------------

ds <- read.csv("Fig7a.csv")
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

#png("Fig7.png",width=6.5,height=4,units="in",res=256)
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
# Fig 8-10
#============================================

ds <- read.csv("Test_n5-100_i0-600mn.csv")
ds

ds$Ihat<-na_if(ds$Ihat,-1)  #may need dplyr

summary(ds[c(4:9)])

summary(ds[c(4:9)])
by(ds[c(4:9)],ds$Group,summary)

#boxplot(ds$Maxdif~ds$Group,varwidth=TRUE,horizontal=TRUE,xlab="dmax",ylab="Run")
#boxplot(ds$Dist~ds$Group,varwidth=TRUE,horizontal=TRUE,xlab="Distance",ylab="Run")
#boxplot(ds$Span~ds$Group,varwidth=TRUE,horizontal=TRUE,xlab="Span",ylab="Run")
#boxplot(ds$IQR~ds$Group,varwidth=TRUE,horizontal=TRUE,xlab="IQR",ylab="Run")
#boxplot(ds$Ihat~ds$Group,varwidth=TRUE,horizontal=TRUE,xlab="Ihat",ylab="Run")
#boxplot(ds$Mid_Date~ds$Group,varwidth=TRUE,horizontal=TRUE,xlab="Year",ylab="Run")

png("Fig8.png",width=6.5,height=4,units="in",res=256)
g1 <- ggplot(ds,aes(y=Maxdif, x=Group)) + labs(x="Run", y="dmax") +
  scale_y_continuous(limits=c(0,600), breaks=seq(0,600,50)) + theme_bw() 
  g1 + geom_boxplot() + coord_flip()
dev.off()
png("Fig8b.png",width=6.5,height=4,units="in",res=256)
g2 <- ggplot(ds,aes(y=Dist, x=Group, )) + labs(x="Run", y="Dist") +
  scale_y_continuous(limits=c(0,600), breaks=seq(0,600,50)) + theme_bw() 
  g2 + geom_boxplot() + coord_flip()
dev.off()
png("Fig8c.png",width=6.5,height=4,units="in",res=256)
g3 <- ggplot(ds,aes(y=Span, x=Group, )) + labs(x="Run", y="Span") +
  scale_y_continuous(limits=c(0,600), breaks=seq(0,600,50)) + theme_bw() 
  g3 + geom_boxplot() + coord_flip()
dev.off()
png("Fig8d.png",width=6.5,height=4,units="in",res=256)
g4 <- ggplot(ds,aes(y=IQR, x=Group, )) + labs(x="Run", y="IQR") +
  scale_y_continuous(limits=c(0,600), breaks=seq(0,6600,50)) + theme_bw() 
  g4 + geom_boxplot() + coord_flip()
dev.off()
png("Fig9.png",width=6.5,height=4,units="in",res=256)
g5 <- ggplot(ds,aes(y=Ihat, x=Group, )) + labs(x="Run", y="Ihat") +
  scale_y_continuous(limits=c(0,600), breaks=seq(0,600,50)) + theme_bw() 
  g5 + geom_boxplot() + coord_flip()
  dev.off()
png("Fig10.png",width=6.5,height=4,units="in",res=256)
g6 <- ggplot(ds,aes(y=Mid_Date, x=Group)) + labs(x="Run", y="Mean") +
  scale_y_continuous(limits=c(975,1225), name="Mean Date", breaks=seq(975,1225,25)) + theme_bw() 
  g6 + geom_boxplot() + coord_flip()
dev.off()

  ds50 <-subset(ds,SD == 50) 
  g1 <- ggplot(ds50,aes(y=Maxdif, x=Group)) + labs(x="Run", y="dmax") +
    scale_y_continuous(limits=c(0,600), breaks=seq(0,600,50)) + theme_bw() 
  g1 + geom_boxplot() + coord_flip()
  g2 <- ggplot(ds50,aes(y=Dist, x=Group, )) + labs(x="Run", y="Dist") +
    scale_y_continuous(limits=c(0,600), breaks=seq(0,600,50)) + theme_bw() 
  g2 + geom_boxplot() + coord_flip()
  g3 <- ggplot(ds50,aes(y=Span, x=Group, )) + labs(x="Run", y="Span") +
    scale_y_continuous(limits=c(0,600), breaks=seq(0,600,50)) + theme_bw() 
  g3 + geom_boxplot() + coord_flip()
  g4 <- ggplot(ds50,aes(y=IQR, x=Group, )) + labs(x="Run", y="IQR") +
    scale_y_continuous(limits=c(0,600), breaks=seq(0,6600,50)) + theme_bw() 
  g4 + geom_boxplot() + coord_flip()
  g5 <- ggplot(ds50,aes(y=Ihat, x=Group, )) + labs(x="Run", y="Ihat") +
    scale_y_continuous(limits=c(0,600), breaks=seq(0,600,50)) + theme_bw() 
  g5 + geom_boxplot() + coord_flip()
  g6 <- ggplot(ds50,aes(y=Mid_Date, x=Group)) + labs(x="Run", y="Mean") +
    scale_y_continuous(limits=c(975,1225), name="Mean Date", breaks=seq(975,1225,25)) + theme_bw() 
  g6 + geom_boxplot() + coord_flip()
  
#  g11 <- ggplot(ds,aes(y=Maxdif, x=Group)) + labs(x="Run", y="dmax") +
#    scale_y_continuous(limits=c(0,500), breaks=seq(0,500,50)) + theme_bw() 
#  g11 + geom_violin() + coord_flip()
    
#  g16<- ggplot(ds,aes(y=Ihat, x=Group, )) + labs(x="Run", y="Ihat") +
#    scale_y_continuous(limits=c(0,500), breaks=seq(0,500,50)) + theme_bw() 
#  g16 + geom_violin() + coord_flip()

  
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
# Fig 11 & 12
#=======================================
copan <-read.csv("Fig11-CopanR.csv")
copan <-data.frame(copan)
summary(copan)
png("Fig11.png",width=6.5,height=4,units="in",res=256)
ggplot(copan, aes(Date)) + #  ggtitle(label) +
  geom_histogram(binwidth=50,breaks=seq(400,1250,50),) + theme_bw() +
#  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_x_continuous(breaks=seq(400,1250,50),limits=c(400,1250),aes(Name="Year")) + 
  scale_y_continuous(breaks=seq(0,350,50,),name="Count") 
dev.off()

C10L2TR <-read.csv("Fig12-10L2TR.csv")
C10L2TR <-data.frame(C10L2TR)
summary(C10L2TR)
png("Fig12.png",width=6.5,height=4,units="in",res=256)
ggplot(C10L2TR, aes(Date)) + #  ggtitle(label) +
  geom_histogram(breaks=seq(550,950,25)) + theme_bw() +
  #  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_x_continuous(breaks=seq(550,950,50),limits=c(550,950),aes(Name="Year")) + 
  scale_y_continuous(breaks=seq(0,12,2),limits=c(0,12),name="Count") 
dev.off()


GP <-read.csv("GP36.csv")
GP <-data.frame(GP)
summary(GP)

#png("Fig13.png",width=6.5,height=4,units="in",res=256)
ggplot(GP, aes(Date)) + #  ggtitle(label) +
  geom_histogram(breaks=seq(300,750,50)) + theme_bw() +
  #  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_x_continuous(breaks=seq(300,750,50),limits=c(300,750),aes(Name="Year")) + 
  scale_y_continuous(breaks=seq(0,12,2),limits=c(0,12),name="Count") 
#dev.off()


GP <-read.csv("GPRP.csv")
GP <-data.frame(GP)
summary(GP)

#png("Fig14.png",width=6.5,height=4,units="in",res=256)
ggplot(GP, aes(Date)) + #  ggtitle(label) +
  geom_histogram(breaks=seq(300,750,50)) + theme_bw() +
  #  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_x_continuous(breaks=seq(300,750,50),limits=c(300,750),aes(Name="Year")) + 
  scale_y_continuous(breaks=seq(0,12,2),limits=c(0,12),name="Count") 
#dev.off()

#=======================================
#Initial Computatins with median, compare with mean
#=======================================
  mn <- read.csv("mean100.csv")
  md <- read.csv("median100.csv")
  
  summary(mn[c(4:9)])
  summary(md[c(4:9)]) 
  
  summary(mn[c(4:9)])
  by(mn[c(4:9)],ds$Group,summary)
  summary(mnd[c(4:9)])
  by(md[c(4:9)],ds$Group,summary)
  
plot(mn$Maxdif,md$Maxdif)

#main difference is in maxdif
#more low outliers with median.
#mean works somewhat better
  
#=======================================
#Euclidean Distance Vs Average Deviation
#=======================================
ed <-  read.csv("Test_ed.csv")
ad <- read.csv("Test_ad.csv")

summary(ed[c(4:9)])
summary(ad[c(4:9)])
  
summary(ed[c(4:9)])
  by(ed[c(4:9)],ed$Group,summary)
summary(ad[c(4:9)])
  by(ed[c(4:9)],ed$Group,summary)
  
# ED works better based on small number of runs, but neither is good 

  
#=======================================
# Cowgill Figure - Fig 1
#=======================================

seedval <- round(runif(1,min=0,max=100000))
#seedval <-  36490  #this one works but remove to retry
set.seed(seedval)
seedval
  
datevec <-  matrix(round(rnorm(40, mean=1000, sd=50)),ncol=1)   
dates <- cbind(1:40,datevec,50)
colnames(dates) <- c("Sample","Year","SD")
dates

datesdf<-data.frame(dates)
datesdf

png("Fig1.png",width=6.5,height=4,units="in",res=256)
label <-"40 Dates Pertaining to the Year 1000"

ggplot(datesdf, aes(Sample, Year)) + theme_bw() +   
  geom_errorbar(aes(ymin=Year-SD, ymax=Year+SD, width=0.25))+
  geom_point(shape=19,size=2)+
  geom_hline(yintercept=1000, size=.8) +
  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_y_continuous(limits=c(800,1200),breaks=seq(800,1200,50)) + 
  coord_flip() + ggtitle(label)

dev.off()

#=======================================
# Cowgill Figure - Fig 2
#=======================================

seedval <- round(runif(1,min=0,max=100000))
#seedval <-  4087  #this one works but remove to retry
set.seed(seedval)
seedval

truedate <- matrix(round(runif(30,min=950,max=1050)))
error <- matrix(round(rnorm(30,mean=0,sd=50)))
measdate <- truedate+error

dates <- cbind(1:30,truedate,error,measdate,50)
colnames(dates) <- c("Sample","True","Error","Measured","SD")
dates

datesdf<-data.frame(dates)
datesdf


png("Fig2.png",width=6.5,height=4,units="in",res=256)
label <-"30 True Dates from the Interval 950-1000 with their Measured Dates"

ggplot(datesdf, aes(Sample, Measured)) +  theme_bw() +
  geom_errorbar(aes(ymin=Measured-SD, ymax=Measured+SD, width=0.25))+
  geom_point(shape=19,size=2)+
  geom_point(aes(y=True,x=Sample), shape=9, size=2) +
  geom_hline(yintercept=950, size=.8) +
  geom_hline(yintercept=1050, size=.8) +
  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_y_continuous(limits=c(800,1200),breaks=seq(800,1200,50)) + 
  coord_flip() + ggtitle(label)

dev.off()

#=======================================
# Cowgill Figure - Fig 3
#=======================================
few <- 5
imin <- 900
imax <- 1100

seedval <- round(runif(1,min=0,max=100000))
#seedval <-  67662  #this one works but remove to retry
set.seed(seedval)
seedval

truedate3 <- matrix(round(runif(few,min=imin,max=imax)))
error3 <- matrix(round(rnorm(few,mean=0,sd=50)))
measdate3 <- truedate3+error3

dates3 <- cbind(1:few,truedate3,error3,measdate3,50)
colnames(dates3) <- c("Sample","True","Error","Measured","SD")
dates3

dates3df<-data.frame(dates3)
dates3df

#note Run this repeatedly until mesured dates are in interval and don't span the whole interval
png("Fig3.png",width=6.5,height=4,units="in",res=256)
label <- paste(few," True Dates from the Interval ",imin,"-",imax," with their Measured Dates",sep="")
subtitle <- "Measured Dates do not Represent the Full range of the Interval"

ggplot(dates3df, aes(Sample, Measured)) +  theme_bw() + 
  geom_errorbar(aes(ymin=Measured-SD, ymax=Measured+SD, width=0.25))+
  geom_point(shape=19,size=2)+
  geom_point(aes(y=True,x=Sample), shape=9, size=2) +
  geom_hline(yintercept=imin, size=.8) +
  geom_hline(yintercept=imax, size=.8) +
  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_y_continuous(limits=c(800,1200),breaks=seq(800,1200,50)) + 
  coord_flip() + ggtitle(label)

dev.off()

#=======================================
# Cowgill Figure - Fig 4 - Normal
#=======================================

#n4 <- ggplot(data.frame(x = c(-4,4)), aes(x = x)) +
#  stat_function(fun = dnorm)
#n4

x <- matrix(c(800:1200),ncol=1)
y <- matrix(dnorm(x,mean=1000,sd=50))

matxy <- cbind(x,y)
colnames(matxy) <- c("Year","P")
matxy <- data.frame(matxy)
matxy

png("Fig4.png",width=6.5,height=3,units="in",res=256)
label <- "Normal Probability Distribution"

ggplot(matxy, aes(Year,P)) + theme_bw() +  
  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_x_continuous(limits=c(800,1200),breaks=seq(800,1200,50)) + 
  scale_y_continuous(breaks=0, minor_breaks=NULL) +
  geom_line() + ggtitle(label)

dev.off()

#=======================================
# Cowgill Figure - Fig 6 - truncated normal
#=======================================

x <- matrix(c(800:900,900:1100,1100:1200),ncol=1)
y <- matrix(dnorm(x,mean=1000,sd=50))

matxy <- cbind(x,y)
colnames(matxy) <- c("Year","P")
matxy <- data.frame(matxy)
matxy

for (i in 1:403) {
  if ((matxy$Year[i]<900) || (matxy$Year[i]>1100)) matxy$P[i] <- 0 
}
matxy[101,2]<-0
matxy[303,2]<-0
matxy

png("Fig6.png",width=6.5,height=3,units="in",res=256)
label <- "Truncated Normal Probability Distribution at +/- 2SD"

ggplot(matxy, aes(Year,P)) + theme_bw() +  
  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_x_continuous(limits=c(800,1200),breaks=seq(800,1200,50)) + 
  scale_y_continuous(breaks=0, minor_breaks=NULL) +
  geom_line() + ggtitle(label)

dev.off()

#=======================================
# Cowgill Figure - Fig 5 - rectangular
#=======================================

x <- rbind(c(800,0),c(950,0),c(950,1),c(1050,1),c(1050,0),c(1200,0))
colnames(x) <- c("Year","P")
x <- data.frame(x)
x

png("Fig5.png",width=6.5,height=2.5,units="in",res=256)
label <- "Rectangular Probability Distribution"

ggplot(x, aes(Year,P)) + theme_bw() + 
  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_x_continuous(limits=c(800,1200),breaks=seq(800,1200,50)) + 
  scale_y_continuous(breaks=0, minor_breaks=NULL, limits=c(0,1.15)) +
    geom_line() + ggtitle(label)

dev.off()

#=======================================
# Cowgill Figure - Fig 7
#=======================================

ndates <- 500
sO <- 100
sM <- 50
sT <- (sO*sO -sM*sM)^0.5
sT

#mode=N2 normal R=rectangular
#mode <- "N2"
mode="R"

if (mode == "N2"){
  K <- 4.55
} else {
  K <- 3.46
}

Ihat <- K*sT
Ihat
sample_sd <- 0

while (sample_sd != 100) {
#seedval <- 314159
#seedval <- round(runif(1,min=0,max=100000))
seedval <-  63617  #this one works but remove to retry
set.seed(seedval)

truedate <- matrix(1:ndates)
error <- matrix(round(rnorm(ndates,mean=0,sd=sM)))

i <- 1

while (i < (ndates+1)) {
 if (mode=="N2") {
    truedate[i] <- round(rnorm(1,mean=1000,sd=sT)) 
    if ((truedate[i]>=(1000-round(Ihat/2))) && (truedate[i]<=(1000+round(Ihat/2)))) i <- i+1  
  } else {
    truedate[i] <- round(runif(1,min=round(1000-Ihat/2),max=round(1000+Ihat/2)))
    print(c(i,truedate[i]))
    i <- i+1
  }
} 

measdate <- truedate+error
sample_sd <- round(sd(measdate))
print(c(seedval,sample_sd))
}

dates <- cbind(1:ndates,truedate,error,measdate)
colnames(dates) <- c("Sample","True","Error","Measured")
datesdf<-data.frame(dates)
datesdf

min(datesdf$Measured)
max(datesdf$Measured)


ggplot(datesdf, aes(Measured)) +   ggtitle(label) +
  geom_histogram(binwidth=50) + theme_bw() +
  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_x_continuous(breaks=seq(700,1300,50)) + 
  scale_y_continuous(name="Count", limits=c(-25,100), breaks=NULL,) 

ggplot(datesdf, aes(True)) +   ggtitle(label) +
  geom_histogram(binwidth=50) + theme_bw() +
  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_x_continuous(breaks=seq(700,1300,50),limits=c(700,1300))+ 
  scale_y_continuous(name="Count", limits=c(-25,100), breaks=NULL) 

# ggplot(datesdf, aes(Sample, Measured)) +  theme_bw() + 
#   geom_errorbar(aes(ymin=Measured-50, ymax=Measured+50, width=0.25))+
#   geom_point(shape=19,size=2)+
#   geom_point(aes(y=True,x=Sample), shape=9, size=2) +
#   geom_hline(yintercept=1000-round(Ihat/2), size=.8) +
#   geom_hline(yintercept=1000+round(Ihat/2), size=.8) +
#   theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
#   scale_y_continuous(limits=c(800,1200),breaks=seq(800,1200,50)) + 
#   coord_flip() + ggtitle(label)


sdm <- rbind(c(900,-5),c(950,-5))
colnames(sdm) <- c("Year","Zero")
sdm <- data.frame(sdm)
sdm

sdo <- rbind(c(1000,-5),c(1100,-5))
colnames(sdo) <- c("Year","Zero")
sdo <- data.frame(sdo)
sdo

rinterval <- rbind(c(1000-300/2,-15),c(1000+300/2,-15))
colnames(rinterval) <- c("Year","Zero")
rinterval <- data.frame(rinterval)
rinterval

ninterval <- rbind(c(1000-403/2,-26),c(1000+403/2,-26))
colnames(ninterval) <- c("Year","Zero")
ninterval <- data.frame(ninterval)
ninterval

lab1  <- expression(sigma[m] == 50)
lab1a <- expression(S[o] == 100)
lab2 <- "Estimated True Interval (Rectangular Distribution)"
lab3 <- "Estimated True Interval (Truncated Normal Distribution)"

png("Fig7.png",width=6.5,height=4,units="in",res=256)
label <- paste(ndates," Measured Dates, So=100 Measurement Sigma=50",sep="")

ggplot(datesdf, aes(Measured)) +   ggtitle(label) +
  geom_histogram(binwidth=50) + theme_bw() +
  theme(plot.title=element_text(size=12), plot.subtitle=element_text(size=12)) +
  scale_x_continuous(breaks=seq(700,1300,50)) + 
  scale_y_continuous(name="Count", limits=c(-32,100), breaks=NULL,) +
  geom_line(data=sdm, aes(x=Year, y=Zero), size=2) +
  geom_text(aes(x = 925, y = -9, label= as.character(lab1)), parse=TRUE,  size=3.5) +
  geom_line(data=sdo, aes(x=Year, y=Zero), size=2) +
  geom_text(aes(x = 1050, y = -9, label= as.character(lab1a)), parse=TRUE, size=3.5) +
  geom_line(data=rinterval, aes(x=Year, y=Zero), size=2) +
  geom_text(aes(x = 1000, y = -19, label=lab2),   size=3.5) +
  geom_line(data=ninterval, aes(x=Year, y=Zero), size=2) +
  geom_text(aes(x = 1000, y = -31, label=lab3),   size=3.5)

dev.off()

ggplot(datesdf, aes(True)) + theme_bw() +
  geom_histogram(binwidth=50)

ggplot(datesdf, aes(Measured)) + 
  geom_dotplot(binwidth=20, dotsize=.5)

ggplot(datesdf, aes(True)) + 
  geom_dotplot(binwidth=20, dotsize=.5)

