\name{entropy}
\alias{entropy}
\title{Entropy}

\description{
Calculate the entropy of data for a given mixture template.}

\usage{
entropy(dat, distr, class, obj)

}
\arguments{

\item{dat}{A data frame.}
\item{distr}{A three letter string indicating component distribution, "mvn"=normal distribution, "mvt"=t-distribution,"msn"=skew normal distribution, "mst"=skew t-ditribution. }
\item{class}{A string, class="jcm", or not.}
\item{obj}{A list of template parameters.}

}
\details{
}
\value{A positive value.}
\references{
M.N., Do and M., Vetterli(2002). Wavelet-Based Texture Retrieval Using Generalized Gaussian 
Density and Kullback-Leibler Distance. IEEE Transactions on Image Processing, 11(2), 146-158.


T.M., Cover and J.A., Thomas(1991). Elements of Information Theory. New York: Wiley.
}
\author{}
\note{
}

\seealso{
\code{\link{do.kld}}
}
\examples{


}


