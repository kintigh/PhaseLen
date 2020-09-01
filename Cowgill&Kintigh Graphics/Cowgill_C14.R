
#One time: Be sure to install packages
#          Tools>Install Packages then add ggplot2,plyr,gmodels,deducer
#
#Each Time:In Packages tab at lower right, check plyr, ggplot2, Psych, labdsv, ca

library("ggplot2")
library("plyr")

#
setwd()   
setwd("c:/users/kintigh/dropbox (ASU)/Archy Research/C14/AmAntiq/Cowgll Figs")

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

#=======================================
