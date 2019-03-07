\name{plotmixt.skew}
\alias{plotmixt.skew}
\title{Overlayed contours}
\description{
It is a function to plot overlayed multiple distribution templates.
}
\usage{
plotmixt.skew(retfiles,which.dim,distr="mvt", 
cont=c(25,50,75), xlim, ylim, zlim, gridsize, 
nrand=1e5, var.label, line.col, marginal.col,  
feature.col, lwd=1,...)
}
\arguments{
\item{retfiles}{The return of dir()}
\item{which.dim}{which dimensions are to be plotted.}
\item{distr}{ A three letter string indicating the type of distribution to be fit. See Details.}
\item{cont}{The levels of contours in percentage.}
\item{xlim}{The range of the first variable}
\item{ylim}{The range of the second variable}
\item{zlim}{The range of the third variable}
\item{gridsize}{The grid size}
\item{nrand}{The number of random number to be generated.}
\item{var.label}{The labels of variables}
\item{line.col}{The colors of contour lines}
\item{marginal.col}{The color of lines}
\item{feature.col}{The color of features}
\item{lwd}{The width of lines}
\item{...}{others}
}
\details{
The source code is stored in R program "JCM_OverlayPlots.R".
}
\value{
}

\seealso{\code{\link{emmix}}
}
\examples{

\dontrun{

setwd("../Desktop/JCM_final")
require(emmix)

splom.settings <- list(layout.heights=list(top.padding=0.5, 
main.key.padding=0, key.axis.padding=0, axis.xlab.padding=0, 
xlab.key.padding=0, key.sub.padding=0, bottom.padding=0, strip=0), 
layout.widths=list(left.padding=0, key.ylab.padding=0,ylab.axis.padding=0, 
axis.key.padding=0,right.padding=0.5))

labels = c("SLP76"  ,"ZAP70",  "CD4",  "CD45RA")

whichcols = c(3,2,4)

retfiles<- paste( c(0,1,3), "min-jca-mvt-5.ret", sep="")

col5=rainbow(5)[c(1,3,5)]

set.seed(8192)

#1. code for "3_timepoints.jpeg"

jpeg(filename="3_timepoints.jpeg", height=512, width=512)

plotmixt.skew(retfiles, which.dim=whichcols, distr="mvt", 
xlim=c(2,10), ylim=c(2,10), zlim=c(2,10), 
par.settings=splom.settings, nrand=1e4, 
var.label=labels[whichcols], lwd=1, cont=c(50,90), 
line.col=col5, add.feature=TRUE, feature.col=col5, 
key=list(space=list("bottom"), 
text=list(c("0 min", "1 min", "3 min")), 
lines=list(lwd=1, col=col5), columns=3))

dev.off()

#2. 3sample plots

dirpath  <- dir("d:/uq/Data/SixPanel",full.names=T)


pid <- 3  #1 min class

files <- dir(dirpath[pid],full.names=T)

dat1 <- read.table(files[1],header=T)
dat2 <- read.table(files[2],header=T)
dat3 <- read.table(files[3],header=T)

dat <- rbind(dat1,dat2,dat3)

p=4
g=5

retAll <- emmix(dat,g,"mvt")

init1 <- emmix(dat1,g,"mvt",init=retAll)
init2 <- emmix(dat2,g,"mvt",init=retAll)
init3 <- emmix(dat3,g,"mvt",init=retAll)


filenames <- dir(dirpath[pid])

dput(init1,paste(filenames[1],".ret",sep=''))
dput(init2,paste(filenames[2],".ret",sep=''))
dput(init3,paste(filenames[3],".ret",sep=''))

dput(retAll,"1min-mvt-5.ret")

#-----------------------------

dirpath  <- dir("d:/uq/Data/SixPanel",full.names=T)

pid <- 3  #1 min class

filenames <- dir(dirpath[pid])

labels = c("SLP76"  ,"ZAP70",  "CD4",  "CD45RA")

whichcols = c(3,2,4)

retfiles<-c(paste(filenames,".ret",sep=''), paste( "1min-mvt-5.ret", sep=""))

retfiles<-c(paste( "1min-mvt-5a.ret", sep=""), paste( "1min-mvt-5.ret", sep=""))

pa <- c(rainbow(5)[c(1,3,5)],"blue")


set.seed(12345)

jpeg("1min(three samples).jpeg",width=512,height=512)

plotmixt.skew(retfiles, whichcols, distr="mvt", 
xlim=c(1.5,8), ylim=c(1.5,8), zlim=c(1.5,10),  
par.settings=splom.settings, nrand=10000,  
var.label=labels[whichcols], lwd=1, 
line.col=pa, cont=c(50,90),  key=list(space=list("bottom"), 
text=list(c("sample1", "sample2", "sample3","template")), 
lines=list(lwd=1, col=pa), columns=4))

dev.off()
}

}
