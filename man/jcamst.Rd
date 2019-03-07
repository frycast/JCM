\name{jcamst}
\alias{jcamst}
\alias{jcamvt}
\title{
Joint clustering and matching
}
\description{
Jointly clustering and matching multiple samples}
\usage{
jcamst(datfiles, p, g, initobj, ncov = 3, dofon = 0, debug = 1, 
itmax = 1000, epsilon = 1e-05, header = TRUE,ftype="csv")
jcamvt(datfiles, p, g, initobj, ncov = 3, dofon = 0, debug = 1, 
itmax = 1000, epsilon = 1e-05, header =  TRUE,ftype="csv")

}
\arguments{
  \item{datfiles}{Full path to where the data are stored.}
  \item{p}{The dimension of data.}
  \item{g}{The number of components.}
  \item{initobj}{The list of initial parametrers.}
  \item{ncov}{The type of covariance structure.}
  \item{dofon}{When set, the degrees of freedom for each component will be same.}
  \item{debug}{Whether to print log-likehood for each iteration.}
  \item{itmax}{The maximum number of iterations.}
  \item{epsilon}{Stopping rule based on relative change in the log-likelihood.}
  \item{header}{Are there headers in the data file.}
\item{ftype}{A string indicating which type of format the data are stored, currently only "txt" and "csv" are handled.}
}
\details{
}
\value{
  \item{error}{Error code, 0 = normal exit;  1 = did not converge within \code{itmax} iterations.}
  \item{loglik}{The log likelihood at convergence.}
  \item{bic}{Bayesian Information Criterion (BIC).}
  \item{pro}{A vector with the mixing proportions.}
  \item{alfa}{A numeric matrix with each column corresponding to the mean.}
  \item{mu}{A numeric matrix with each column corresponding to the mean.}
  \item{sigma}{An array of dimension (p,p,g) with first two arguments corresponding to the covariance matrix of each component.}
  \item{dof}{A vector of degrees of freedom for each component.}
  \item{delta}{A p by g matrix with each column corresponding to a skew parameter vector.}
  \item{theta}{A matrix with each column corresponding to the variances of the random effects iterm a.}
  \item{thetu}{A vector with the variance of the random effects item b.}
  \item{distr}{A three letter string standing for distribution.}
  \item{class}{Class should be "jcm".}

}
\references{
}
\author{
}
\note{
}


\seealso{
\code{\link{msmvt}}
\code{\link{msmvn}}
\code{\link{emmix.flow}}
\code{\link{jcamst.predict}}

}
\examples{

data(Lympho)

# we calculate the template of class '0min'.

datnames <- (Lympho$names)[1:5]

# make a temp dir for data files

system("mkdir templym")

datfiles <- paste("templym/",datnames,".txt",sep='')

for(id in 1:5) {

S   <- (Lympho$data)[[id]]
write.table(S,datfiles[id],row.names=FALSE)

}

# set initial values

init<- Lympho$init

p <- ncol(S)
g <- length(init$pro)

init$theta <- array(0.1,c(p,g))
init$thetu <- rep(0.1,g)

#do individuals
\dontrun{
obj <- jcamvt(datfiles[1],p,g,init,itmax=100,ftype="txt")
obj <- jcamst(datfiles[1],p,g,init,itmax=100,ftype="txt")
}

# do class template
\dontrun{

obj <- jcamvt(datfiles,p,g,init,itmax=100,ftype="txt")
obj <- jcamst(datfiles,p,g,init,itmax=100,ftype="txt")

}

# Do JCM MST alternatively
\dontrun{


id=1

S   <- (Lympho$data)[[id]]

p=4
g=4
n <- nrow(S)

sl <- sample(1:n,min(2000,n))
x11()
emmix.flow(S[sl,],diag.panel=TRUE)

init <- emmix(S[sl,],g,"mst",itmax=100)

init$distr="mst"

emmix.flow(S[sl,],init,diag.panel=TRUE)

init$theta<- array(0.1,c(p,g))
init$thetu<- rep(0.1,g)

obj <- jcamst(datfiles[1],p,g,init,itmax=100,ftype="txt")

pre <- jcamst.predict(S, g, obj$pro, obj$mu, obj$sigma, 
obj$dof, obj$delta,obj$theta, obj$thetu)

#plot the aligned data
x11()
emmix.flow(pre$eee,diag.panel=TRUE)
title(main = list(paste("Aligned flow sample:",datnames[id]), 
cex=1.5,col="blue", font=3))

#plot the data
x11()
emmix.flow(S,clust=pre$clust,diag.panel=TRUE)
title(main = list(paste("flow sample:",datnames[id]), 
cex=1.5,col="blue", font=3))

}

}

