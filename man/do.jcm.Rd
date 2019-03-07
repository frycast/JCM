\name{do.jcm}
\alias{do.jcm}
\title{Fit JCM on Each Sample of Panel}
\description{
It is a pipeline function to fit each sample of a panel into JCM model.
}
\usage{
do.jcm(data="4min",panel="panel1",samples,p=4,g=2,distr="mvt",itmax=100,
header=TRUE,ftype="CSV")
}
\arguments{
\item{data}{A string of the data name, which refers to a location where the panel data are stored.}
\item{panel}{A string of the panel name, which refers to a subfolder where the real data are stored.}
\item{samples}{The vector of sample names}
\item{p}{The number of variables of data}
\item{g}{The number of components of the mixture model}
\item{distr}{ A three letter string indicating the type of distribution to be fit. See Details.}
\item{itmax}{The maximum of iterations}
\item{header}{A logical value, indicating whether there is a header in the data file.}
\item{ftype}{A string indicating which type of format the data are stored, currently only "txt" and "csv" are handled.}
}
\details{
This function has to be run after domix,because the domix results are used as inital values.
}
\value{
A ret file is returned for each sample. All estimated information of a sample are stored in its own ret file and can be accessed via dget() function in R. 
}

\seealso{\code{\link{do.mix}},\code{\link{do.template}},\code{\link{do.unzip}}
}
\examples{

\dontrun{

samples <- c(info[info[,2]=="LNP_lo",1])

do.jcm(data="4min",panel="panel1",samples,p=4,g=2,
distr="mvt",header=TRUE,ftype="CSV")

}

}
