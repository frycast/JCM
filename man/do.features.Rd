\name{do.features}
\alias{do.features}
\title{Extract the feature matrix}
\description{
It is a pipeline function to get the features for each individual samples.
}
\usage{
do.features(retpath='.',p=4,g=2,distr="mst",
clusters=c("base","mould"),data="0min",
panel="panel1",group="LNP_hi",samples,
markers=c("p-SFK","BCL2" ,"CD20","p-ERK"),class="jcm") 

}
\arguments{
\item{retpath}{A string of path to where the return files are stored.}
\item{p}{The number of variables of data}
\item{g}{The number of components of the mixture model}
\item{distr}{ A three letter string indicating the type of distribution to be fit. See Details.}
\item{clusters}{A vector of strings for names of the clusters.}
\item{data}{A string of the data name, which refers to a location where the panel data are stored.}
\item{panel}{A string of the panel name, which refers to a subfolder where the real data are stored.}
\item{group}{A string of clinical group.}
\item{samples}{A vector of sample names.}
\item{markers}{A vector of marker names.}
\item{class}{A string of class, either "jcm" or space.}
}
\details{

}
\value{
"batch_features_data_panel_group.txt"
"batch_features_data_panel_group.gct"

}

\seealso{\code{\link{do.jcm}},\code{\link{do.mix}},\code{\link{do.unzip}}
}
\examples{

\dontrun{

info <- as.matrix(read.csv("JCM_clinical.csv"))

samples <- c(info[info[,2]=="sensitive",1])

data   <- "0min"

panel  <- "panel1"

group  <- "LNP_lo"

markers<- c("p-SFK","BCL2" ,"CD20","p-ERK")

do.features(retpath='.',p=4,g=2,distr="mst",
clusters=c("base","mould"),
data,panel,group,samples,markers,class="jcm") 

}
}

