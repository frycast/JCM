\name{do.mould}
\alias{do.mould}
\title{Extract the mould component}
\description{
It is a pipeline function to get the mould component for each individual samples.
}
\usage{
do.mould(data,panel,sample,p=4,g=2,distr="mst",class="") 

}
\arguments{
\item{data}{A string of the data name, which refers to a location where the panel data are stored.}
\item{panel}{A string of the panel name, which refers to a subfolder where the real data are stored.}
\item{sample}{A string of sample name.}
\item{p}{The number of variables of data}
\item{g}{The number of components of the mixture model}
\item{distr}{ A three letter string indicating the type of distribution to be fit. See Details.}
\item{class}{A string of class, either "jcm" or space.}
}
\details{
This example is to illustrate how to plot the xyplot for multiple samples.
}
\value{
A matrix of mould data.

}

\seealso{\code{\link{do.jcm}},\code{\link{do.mix}},\code{\link{do.unzip}}
}
\examples{

\dontrun{

info <- as.matrix(read.csv("JCM_clinical.csv"))

samples  <- c(info[info[,2]=="sensitive",1])

samples


# [1] "lp-j105"  "lp-j106"  "lp-j108"  "lp-j111"  "lp-j112"  "lp-j116" 
# [7] "lp-j117"  "lp-j120"  "lp-j122"  "lp-j123"  "lp-j124"  "lp-j126" 
#[13] "lp-j127"  "lp-j131a" "lp-j133a" "lp-j134"  "lp-j140"  "lp-j141" 
# these are LNP- samples



data<- "4min"

panel<- "panel4"

distr<-"mst"

class<-""


#Row 1. 105, 108, 111, 112, 116;

lab3<- c("lp-j105","lp-j108","lp-j111","lp-j112","lp-j116")

z1=do.mould(data,panel,"lp-j105",p=4,g=2,distr,class)
z2=do.mould(data,panel,"lp-j108",p=4,g=2,distr,class)
z3=do.mould(data,panel,"lp-j111",p=4,g=2,distr,class)
z4=do.mould(data,panel,"lp-j112",p=4,g=2,distr,class)
z5=do.mould(data,panel,"lp-j116",p=4,g=2,distr,class)



#Row 2. 124, 126, 127, 134, 140;

lab2<- c("lp-j124","lp-j126","lp-j127","lp-j134","lp-j140")

y1=do.mould(data,panel,"lp-j124",p=4,g=2,distr,class)
y2=do.mould(data,panel,"lp-j126",p=4,g=2,distr,class)
y3=do.mould(data,panel,"lp-j127",p=4,g=2,distr,class)
y4=do.mould(data,panel,"lp-j134",p=4,g=2,distr,class)
y5=do.mould(data,panel,"lp-j140",p=4,g=2,distr,class)

samples  <- c(info[info[,2]=="insensitive",1])

samples

# [1] "lp-j101"  "lp-j102"  "lp-j103"  "lp-j109"  "lp-j110"  
#"lp-j114" [7] "lp-j119"  "lp-j121"  "lp-j125"  "lp-j128a"
# these are LNP+ samples


#Row 3. 109, 119, 128a, 121, 103

x1=do.mould(data,panel,"lp-j103",p=4,g=2,distr,class)
x2=do.mould(data,panel,"lp-j109",p=4,g=2,distr,class)
x3=do.mould(data,panel,"lp-j119",p=4,g=2,distr,class)
x4=do.mould(data,panel,"lp-j121",p=4,g=2,distr,class)
x5=do.mould(data,panel,"lp-j128a",p=4,g=2,distr,class)

lab1=c("lp-j103","lp-j109","lp-j119","lp-j121","lp-j128a")

dat <- rbind(x1 ,x2 ,x3 ,x4 ,x5 ,
y1 ,y2 ,y3 ,y4 ,y5 ,z1 ,z2 ,z3 ,z4 ,z5 )

markers <- c("p.PLCg2",  "BCL2", "CD20", "p.STAT5")

colnames(dat) <-markers

labels <- factor(rep(1:15,
c(nrow(x1 ),nrow(x2 ),nrow(x3 ),nrow(x4 ),nrow(x5 ),
  nrow(y1 ),nrow(y2 ),nrow(y3 ),nrow(y4 ),nrow(y5 ),
  nrow(z1 ),nrow(z2 ),nrow(z3 ),nrow(z4 ),nrow(z5 ))),
  labels=c(lab1,lab2,lab3))

#-----------


jpeg("xyplotLPN_final(2).jpeg",height=2048*0.35,width=2048*.5)

xyplot(`p.PLCg2`~`CD20`|labels,layout=c(5,3),data=dat,
panel=function(x,y,...) { 
xrange=range(dat[,1])+c(-0.65,0.5)
yrange=range(dat[,2])+c(-0.65,0.5)
Lab.palette <- colorRampPalette(c("blue", "orange", "red"),
space = "Lab")
panel.smoothScatter(x,y, colramp = Lab.palette,
range.x=list(xrange,yrange),bandwidth=c(0.1,0.1),...)})


dev.off()

#-----------------------------------------

}

}