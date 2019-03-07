\name{jcamst.predict}
\alias{jcamst.predict}
\alias{jcamvt.predict}
\title{Calculate the aligned values (via JCAMST)}
\description{
Using the estimated parameters to calculate the aligned values at each data points, as well as the posterior probability and final partitions of the data. 
}
\usage{
jcamst.predict(dat,g,pro,mu,sigma,dof,delta,theta,thetu)
jcamvt.predict(dat,g,pro,mu,sigma,dof,theta,thetu)

}
\arguments{
  \item{dat}{An n by p numeric matrix, the dataset}
  \item{g}{The number of clusters}
  \item{pro}{A vector of mixing proportions, see Details.}
  \item{mu}{A numeric matrix with each column corresponding to the mean, see Details.}
  \item{sigma}{An array of dimension (p,p,g) with first two dimension corresponding covariance matrix of each component, see Details.}
  \item{dof}{A vector of degrees of freedom for each component, see Details.}
  \item{delta}{A p by g matrix with each column corresponding to a skew parameter vector, see Details.}
  \item{theta}{A p by g matrix with each column corresponding to the variance of the sacling ramdom effect.}
  \item{thetu}{A vector of variance for each component.}
}
\value{
\code{eee} is a n by p matrix giving the aligned values; 
\code{tau} is a n by g matrix of the posterior probability; 
\code{clust} is a vector of the clustering.
}
\seealso{
\code{\link{jcamvt}}.
\code{\link{jcamst}}.
}
\examples{
\dontrun{

# get the data

dat <- read.table("mydata.txt",header=T) 

# get the estimated parameters of JCAMVT model
obj<- dget("xxx.ret")

# calculate the predicted values


pred <- predict.jcamst(dat,length(obj$pro),obj$pro,obj$mu,
obj$sigma,obj$dof,obj$delta,obj$theta,obj$thetu)

# plot the clustering results in 2D pairs
# using the aligned values
plot((pred$eee)[,c(1,2)],col=pred$clust)


}

}
\keyword{cluster}
