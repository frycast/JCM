\name{do.mix}
\alias{do.mix}
\title{Fit the Multivariate Mixture Models on Each Sample of Panel}
\description{
It is a pipeline function to fit each sample of a panel into specified mixture model.
}
\usage{
do.mix(data="4min",panel="panel1",samples,p=4,g=2,distr="mvt",
init=NULL,itmax=100,header=TRUE,ftype="CSV")
}
\arguments{
\item{data}{A string of the data name, which refers to a location where the panel data are stored.}
\item{panel}{A string of the panel name, which refers to a subfolder where the real data are stored.}
\item{samples}{A vector of sample names}
\item{p}{The number of variables of data}
\item{g}{The number of components of the mixture model}
\item{distr}{ A three letter string indicating the type of distribution to be fit. See Details.}
\item{itmax}{The maximum of iterations}
\item{init}{The initial values}
\item{header}{A logical value, indicating whether there is a header in the data file.}
\item{ftype}{A string indicating which type of format the data are stored, currently only "txt" and "csv" are handled.}
}
\details{
The distribution type, determined by the \code{distr} parameter, which may take any one of the following values:
"mvn" for a multivariate normal, "mvt" for a multivariate t-distribution, "msn" for a multivariate skew normal distribution and "mst" for a multivariate skew t-distribution.

This function fit the specified mixtrue model to all samples in the panel of data.
}
\value{
A ret file is returned for each sample. All estimated information of a sample are stored in its own ret file and can be accessed via dget() function in R. 
}

\seealso{\code{\link{do.jcm}},\code{\link{do.template}},\code{\link{do.unzip}}
}
\examples{

\dontrun{

samples <- c(info[info[,2]=="LNP_lo",1])

do.mix(data="4min",panel="panel1",samples,p=4,g=2,
distr="mvt",header=TRUE,ftype="CSV")

}
}
