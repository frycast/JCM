
ddmvn<-function(dat,n,p,mean = rep(0,p),cov = diag(p)                   )
{
exp(ddmix(dat,n,p,1,"mvn", mean,cov,0,rep(0,p)))

}

ddmvt<-function(dat,n,p,mean = rep(0,p),cov=diag(p),nu=4                )
{
exp(ddmix(dat,n,p,1, "mvt", mean,cov,nu,rep(0,p)))
}

ddmsn<-function(dat,n,p,mean=rep(0,p),cov=diag(p),       del = rep(0,p))
{
exp(ddmix(dat,n,p,1, "msn", mean,cov,0,del))
}

ddmst<-function(dat,n,p,mean = rep(0,p),cov=diag(p),nu=4,del = rep(0,p))
{
exp(ddmix(dat,n,p,1, "mst", mean,cov,nu,del))
}



ddmix <- function(dat,n,p,g,distr, mu, sigma, dof=NULL, delta=NULL)
{

if(is.null(dof))
dof <- rep(4,g)

if(is.null(delta))
delta <- array(0,c(p,g))

ndist<-switch(tolower(distr),"mvn"=1,"mvt"=2,"msn"=3,"mst"=4,5)

if(ndist>4) 
stop("the model specified is not available yet")

dat<-as.matrix(dat);

if(n == 1 & (ncol(dat) ==1))
dat<-t(dat)

if(nrow(dat)!=n | ncol(dat)!=p )
stop("dat does not match n and p.")

#is mu,sigma,dof,delta specified correctly?

if(length(c(mu)) != (p*g))
stop(paste("mu should be a ",p, 'by', g, "matrix!"))

if(length(c(sigma)) != (p*p*g))
stop(paste("sigma should be a ",p, 'by', p,'by', g, " array!"))

if(length(c(dof)) != g)
stop(paste("dof should be a ",g, " vector!"))

if(length(c(delta)) != (p*g))
stop(paste("delta should be a ",p, 'by', g, " array!"))

obj<-.C('ddmix',PACKAGE="JCM",
as.double(dat),as.integer(n),
as.integer(p),as.integer(g),as.integer(ndist),
as.double(mu),as.double(sigma),
as.double(dof),as.double(delta),
den = double(n*g),error = integer(1))[10:11] 

if(obj$error) stop("error")

(matrix(obj$den,ncol=g))

}

