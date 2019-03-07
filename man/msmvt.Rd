\name{msmvt}
\alias{msmvt}
\alias{msmvn}
\alias{msmsn}
\alias{msmst}
\title{
Distribution template
}
\description{
Calculate the distribution template of a class of multiple samples}
\usage{
msmvt(datfiles, p, g, initobj, ncov = 3, common = 0, debug = 1, 
itmax = 1000, epsilon = 1e-05, header =  TRUE,ftype="txt")

msmvn(datfiles, p, g, initobj, ncov = 3, debug = 1, 
itmax = 1000, epsilon = 1e-05, header = TRUE,ftype="txt")

msmsn(datfiles, p, g, initobj, ncov = 3, debug = 1, 
itmax = 1000, epsilon = 1e-05, header = TRUE,ftype="txt")

msmst(datfiles, p, g, initobj, ncov = 3, debug = 1, 
itmax = 1000, epsilon = 1e-05, header = TRUE,ftype="txt")
}
\arguments{

\item{datfiles}{Full path to where the data are stored.}
  \item{p}{Dimension of data.}
  \item{g}{The number of components.}
  \item{initobj}{The list of initial parametrers.}
  \item{ncov}{The type of covariance structure.}
  \item{common}{When set, the degrees of freedom for all components will be same.}
  \item{debug}{When set, the log likehood for each iteration will be printed.}
  \item{itmax}{The maximum number of iterations.}
  \item{epsilon}{Stopping rule based on the relative change in the log likelihood.}
  \item{header}{Are there headers in the data files.}
\item{ftype}{The type of data file; default is "txt", alternative is "csv"}
}
\details{
}
\value{
  \item{error}{Error code, 0 = normal exit;  1 = did not converge within \code{itmax} iterations.}
  \item{loglik}{The log likelihood at convergence.}
  \item{bic}{Bayesian Information Criterion (BIC).}
  \item{pro}{A vector with the mixing proportions.}
  \item{mu}{A numeric matrix with each column corresponding to the mean.}
  \item{sigma}{An array of dimension (p,p,g) with first two arguments corresponding to the covariance matrix of each component.}
  \item{dof}{A vector of degrees of freedom for each component, see Details.}
  \item{delta}{A p by g matrix of zeros.}
  \item{distr}{A three letter string indicating the distribution.}
}
\references{
McLachlan G.J. and Krishnan T. (2008). The EM Algorithm and Extensions (2nd). Hoboken, New Jersey: Wiley.

McLachlan G.J. and Peel D. (2000). Finite Mixture Models. New York: Wiley.

}
\author{
}
\note{
}


\seealso{
\code{\link{emmix.flow}}
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

init <- Lympho$init

# do class template

\dontrun{

obj1 <- msmvt(datfiles,p,g,init,itmax=100)

obj2 <- msmvn(datfiles,p,g,init,itmax=100)
}

}
