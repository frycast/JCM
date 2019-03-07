\name{error.rate}
\alias{error.rate}
\alias{F.Measures}
\title{Error Rate}
\description{
Calculate the Error Rate of a partition
}
\usage{
error.rate(clust1,clust2)
F.Measures(Clabels,Klabels)
}
\arguments{
\item{clust1}{An integer vector of cluster label 1 }
\item{clust2}{An integer vector of cluster label 2}
\item{Clabels}{An integer vector of the true membership labels}
\item{Klabels}{An integer vector of the predicted labels}

}
\details{
clust1 and clust 2 must match, i.e, same number of clusters
}
\value{
error.rate gives Error Rate
}
\examples{
clu1<-c(1,2,3,1,1,2,2,3,3)
clu2<-c(2,2,2,1,1,1,3,3,3)
error.rate(clu1, clu2)
}
\keyword{cluster}
\keyword{datasets}
