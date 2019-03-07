\name{do.kld}
\alias{do.kld}
\title{Approximate Kullback-Leibler Distance}
\description{
Calculate the Approximate Kullback-Leibler Distance (KLD).
}
\usage{
do.kld(datfile, SelfTemplate, ClassTemplates, st, 
distr = "mvt", class = "jcm",header=TRUE,ftype="csv")
}
\arguments{

\item{datfile}{A string where the data are stored.}
\item{SelfTemplate}{A string where the individual template parameters are stored.}
\item{ClassTemplates}{A string where the calss template parameters are stored.}
\item{st}{A vector for the indices of the variables.}
\item{distr}{A three letter string indicating component distribution: "mvn"=normal distribution, "mvt"=t-distribution,"msn"=skew normal distribution, "mst"=skew t-ditribution. }
\item{class}{A string, class="jcm", or not.}
\item{header}{Are there headers in the data file.}
\item{ftype}{The type of data file; default is "csv", alternative is "txt"}
}
\details{
}
\value{
A vector of KLD values
}
\references{
M.N., Do and M., Vetterli(2002). Wavelet-Based Texture Retrieval 
Using Generalized Gaussian Density and Kullback-Leibler Distance. 
IEEE Transactions on Image Processing, 11(2), 146-158.


T.M., Cover and J.A., Thomas(1991). Elements of Information Theory. New York: Wiley.
}
\author{}
\note{}

\seealso{
\code{\link{entropy}}
}

\examples{

data(Lympho)

id=1

datname <- (Lympho$names)[id]

S   <- (Lympho$data)[[id]]

obj <- (Lympho$mvt)[[id]]

#
datfile <- paste(datname,".txt",sep='')
write.table(S,datfile,row.names=FALSE)

retfile <- paste(datname,".ret",sep='')
dput(obj,retfile)


tmpfile <- paste((Lympho$template$names)[1:2],".ret",sep='')

for(j in 1:2) {
obj <- (Lympho$template$template)[[j]]
dput(obj,tmpfile[j])
}


st <- 1:4

do.kld(datfile,retfile,tmpfile,st, distr = "mvt", class = "",ftype="txt")

}

