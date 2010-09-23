#!/usr/bin/Rscript
# R script to plot sliding windows of various statistics along a genome
# to be used in conjunction with output from rhomulus.pl
# @author Aaron Darling
# @license GPL
#

library("zoo")

bbb <- read.table("rho_summary.txt")
ccc <- read.table("theta_summary.txt")
ddd <- read.table("delta_summary.txt")
eee <- read.table("tmrca_summary.txt")
fff <- read.table("rhotheta_summary.txt")
ggg <- read.table("rhopersite_summary.txt")
hhh <- read.table("rhodelta_summary.txt")
iii <- read.table("roverm_summary.txt")

meanie <- function(x){ sum(x[!is.na(x)])/sum(!is.na(x)) }
putter <- function(x,y,z,sitenorm){
	blob <- vector()
	for(i in 1:length(x))
	{
		if(sitenorm==TRUE)
		{
		        blob[seq(x[i],y[i],by=1)] <- z[i]/(y[i]-x[i])
       	 	}else{
			blob[seq(x[i],y[i],by=1)] <- z[i]
		}
	}
	blob
}

polyfillNA <- function(lowline,highline,mycolor)
{
i <- 1
prev <- 1
for( i in 1:length(lowline) )
{
	if(is.na(lowline[i]))
	{
		# time to plot
		if(prev < i)
		{
			ii <- i-1
			polygon(x=c(attr(lowline[prev:ii],'index'),rev(attr(highline[prev:ii],'index'))),y=c(as.numeric(lowline[prev:ii]),rev(as.numeric(highline[prev:ii]))), density=NA, col=mycolor)
		}
		prev <- i+1
	}
}
# get the last one
if(!is.na(prev))
{
	ii <- length(lowline)
	polygon(x=c(attr(lowline[prev:ii],'index'),rev(attr(highline[prev:ii],'index'))),y=c(as.numeric(lowline[prev:ii]),rev(as.numeric(highline[prev:ii]))), density=NA, col=mycolor)
}
}

plotGenomeDistribution <- function(dataset,windowsize,windowstep,sitenorm,ptitle,pylab)
{
meanval <- dataset$V2
leftends <- dataset$V4
rightends <- dataset$V5
minval <- dataset$V7
val10 <- dataset$V8
val25 <- dataset$V9
val50 <- dataset$V10
val75 <- dataset$V11
val90 <- dataset$V12
val100 <- dataset$V13


b0 <- putter(leftends,rightends,minval,sitenorm)
b0r <- rollapply(as.zoo(b0),windowsize,meanie,by=windowstep)

b10 <- putter(leftends,rightends,val10,sitenorm)
b10r <- rollapply(as.zoo(b10),windowsize,meanie,by=windowstep)

b25 <- putter(leftends,rightends,val25,sitenorm)
b25r <- rollapply(as.zoo(b25),windowsize,meanie,by=windowstep)

b50 <- putter(leftends,rightends,val50,sitenorm)
b50r <- rollapply(as.zoo(b50),windowsize,meanie,by=windowstep)

b75 <- putter(leftends,rightends,val75,sitenorm)
b75r <- rollapply(as.zoo(b75),windowsize,meanie,by=windowstep)

b90 <- putter(leftends,rightends,val90,sitenorm)
b90r <- rollapply(as.zoo(b90),windowsize,meanie,by=windowstep)

b100 <- putter(leftends,rightends,val100,sitenorm)
b100r <- rollapply(as.zoo(b100),windowsize,meanie,by=windowstep)

bmean <- putter(leftends,rightends,dataset$V2,sitenorm)



plot(log(b50r),ylim=log(c(min(b0r[!is.na(b0r)]),max(b100r[!is.na(b100r)]))),lwd=0.5,type="l",xlab="E. coli K12 MG1655",main=ptitle,ylab=pylab)
polyfillNA(log(b0r),log(b100r),rgb(.9,.9,.9))
polyfillNA(log(b10r),log(b90r),rgb(.75,.75,.75))
polyfillNA(log(b25r),log(b75r),rgb(.5,.5,.5))
polyfillNA(log(b50r),log(b50r),rgb(0,0,0))

lines(x=c(2750000,3920000), y=log(rep(meanie(bmean[2750000:3920000]),2)),col="orange",lwd=0.5)
lines(x=c(1600000,2750000), y=log(rep(meanie(bmean[1600000:2750000]),2)),col="orange",lwd=0.5)
lines(x=c(450000,1600000), y=log(rep(meanie(bmean[450000:1600000]),2)),col="orange",lwd=0.5)
lines(x=c(3920000,4650000), y=log(rep(meanie(c(bmean[1:450000],bmean[3920000:length(bmean)])),2)),col="orange",lwd=0.5)
lines(x=c(1,450000), y=log(rep(meanie(c(bmean[1:450000],bmean[3920000:length(bmean)])),2)),col="orange",lwd=0.5)

}

pdf("rhopersite_sliding.pdf",width=10,height=4.5)
plotGenomeDistribution(ggg,20000,5000,FALSE,"Smoothed posterior distribution of rho per site","log rho per site")
dev.off()
pdf("rhodelta_sliding.pdf",width=10,height=4.5)
plotGenomeDistribution(hhh,20000,5000,FALSE,"Smoothed posterior distribution of rhopersite*delta","log rpsdelta")
dev.off()
pdf("roverm_sliding.pdf",width=10,height=4.5)
plotGenomeDistribution(iii,20000,5000,FALSE,"Smoothed posterior distribution of r/m","log r/m")
dev.off()
pdf("rhotheta_sliding.pdf",width=10,height=4.5)
plotGenomeDistribution(fff,20000,5000,FALSE,"Smoothed posterior distribution of rho/theta","rho/theta")
dev.off()
pdf("rho_sliding.pdf",width=10,height=4.5)
plotGenomeDistribution(bbb,20000,5000,TRUE,"Smoothed posterior distribution of rho per site","rho per site")
dev.off()
pdf("theta_sliding.pdf",width=10,height=4.5)
plotGenomeDistribution(ccc,20000,5000,TRUE,"Smoothed posterior distribution of theta per site","theta per site")
dev.off()
pdf("delta_sliding.pdf",width=10,height=4.5)
plotGenomeDistribution(ddd,20000,5000,FALSE,"Smoothed posterior distribution of delta","delta")
dev.off()
pdf("tmrca_sliding.pdf",width=10,height=4.5)
plotGenomeDistribution(eee,20000,5000,FALSE,"Smoothed posterior distribution of tmrca","tmrca")
dev.off()

