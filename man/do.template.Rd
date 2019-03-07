\name{do.template}
\alias{do.template}
\alias{JCM.template}
\title{Find the distribution template}
\description{
It is a pipeline function to get the distribution template for a particular group of samples. }
\usage{
do.template(data="4min",panel="panel1",group="LNP_lo",
distr="mvt",p=4,g=6,info,itmax=100,header=TRUE,ftype="csv") 
JCM.template(data="4min",panel="panel1",group="LNP_lo",
distr="mvt",p=4,g=6,info,itmax=100,header=TRUE,ftype="csv") 
}
\arguments{
\item{data}{A string of the data name, which refers to a location where the panel data are stored.}
\item{panel}{A string of the panel name, which refers to a subfolder where the real data are stored.}
\item{group}{A string of clinical group.}
\item{p}{The number of variables of data}
\item{g}{The number of components of the mixture model}
\item{distr}{ A three letter string indicating the type of distribution to be fit. See Details.}
\item{info}{The matrix of a two-column data, with variables as patient and condition group.}
\item{itmax}{The maximum of iterations}
\item{header}{A logical value, indicating whether there is a header in the data file.}
\item{ftype}{A string indicating which type of format the data are stored, currently only "txt" and "csv" are handled.}
}
\details{
All samples are used to get the distribution template, 
while parts of samples are used to get initial values. 
JCM.template has to be run after function do.template 
since JCM.tmeplate input the results of do.template as inital values. 
All tempalte information are stored in the returned ret files 
which can be accessed by dget() in R. 
}
\value{
Two ret files are returned for each clinical group. 
One leads with "init", one with "template". 
}

\seealso{\code{\link{do.jcm}},\code{\link{do.mix}},\code{\link{do.unzip}}
}
\examples{

\dontrun{

info <- read.csv("JCM_clinical.csv")

do.template(data="4min",panel="panel1",group="LNP_lo",
distr="mvt",p=4,g=6,info=info,header=TRUE,ftype="csv") 

JCM.template(data="4min",panel="panel1",group="LNP_lo",
distr="mvt",p=4,g=6,info=info,header=TRUE,ftype="csv") 

}
}

