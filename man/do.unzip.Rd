\name{do.unzip}
\alias{do.unzip}
\title{Unzip the panel data}
\description{
It is a pipeline function to unzip the panel data into the project folder.
}
\usage{
do.unzip(data,path='.',exdir='.') 
}
\arguments{
\item{data}{A string refering to the zipped data name.}
\item{path}{A string refering to a location where the zip files are stored.}
\item{exdir}{A string refering to a where the unzipped data go..}
}
\examples{

\dontrun{

data <- "TCR-6classes"
do.unzip(data,path='.',exdir='.') 

}

}
