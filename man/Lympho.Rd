\name{Lympho}
\alias{Lympho}
\docType{data}
\title{
Lymphocytes
}
\description{
Lymopho consists of ten samples of flow cytometric data, as well as fitting information for msmvt and jcamvt.}
\usage{data(Lympho)}
\format{

A list of 6 components:

$names: names of the ten samples;  

$data: list of ten data frames;   

$mvt:  list of MVT templates for the ten individual samples;

$init: initial values;

$jcm: list of JCM templates for the ten individual samples;

$template: a list with two components, $names and $template; 
$names gives four class template names,"0min-mvt", "5min-mvt", "0min-jcm" and "5min-jcm";
$template is the list of four class templates.

}
\details{
}
\source{
http://www.broadinstitute.org/cancer/software/genepattern/modules/FLAME/published_data.html}
\references{
}
\examples{
data(Lympho)


}
\keyword{datasets}
