
rdmvn<-function (n, p,mean = rep(0,p), cov = diag(p)) 
{
#it is a interface function to call mvtnorm::rmvnorm
cov<-as.matrix(cov)
    
if (nrow(cov) != ncol(cov)) {
        stop("cov must be a square matrix")
    }
    if (length(mean) != nrow(cov)) {
        stop("mean and cov have non-conforming size")
    }

rmvnorm(n, mean = mean, sigma = cov,method="chol")

}


rdmvt<-function(n,p,mean = rep(0,p),cov=diag(p),nu=3)
{
cov<-as.matrix(cov)
u<-rgamma(n,nu/2,nu/2)
t(t(rdmvn(n,p,cov=cov)/sqrt(u))+mean)
}


rdmsn<-function(n,p,mean=rep(0,p),cov=diag(p),del=rep(0,p))
{

x<-rdmvn(n,p,mean,cov)
z<-abs(rnorm(n))
as.matrix(z%*%t(del)+x)

}



rdmst<-function(n,p,mean = rep(0,p), cov=diag(p),nu=10,del = rep(0,p))
{
u<-rgamma(n,nu/2,nu/2)
x<-t(t(rdmvn(n,p,cov=cov)/sqrt(u))+mean)
z<-abs(rnorm(n)/sqrt(u))
as.matrix(z%*%t(del)+x)
}


rdemmix2<-function(n,p,g,distr,pro,mu,sigma,dof=NULL,delta=NULL)
{

n0 <- table(sample(1:g, n, replace = TRUE, prob = pro))
nn <- n0
if(length(nn) <g) {
nn <- rep(0,g)
for(i in as.numeric(names(n0)))
nn[i] <- n0[paste(i)]
}

names(nn) <- NULL

rdemmix(nn,p,g,distr,mu,sigma,dof,delta)

}

rdemmix3<-function(n,p,g,distr,pro,mu,sigma,dof=NULL,delta=NULL)
{

if(length(pro) != g)
stop(paste("pro should be a ",g, " vector!"))

n0 <- table(sample(1:g, n, replace = TRUE, prob = pro))
nn <- n0
if(length(nn) <g) {
nn <- rep(0,g)
for(i in as.numeric(names(n0)))
nn[i] <- n0[paste(i)]
}

names(nn) <- NULL

dat <- rdemmix(nn,p,g,distr,mu,sigma,dof,delta)

list(data = dat, cluster = rep(1:g,nn) )

}


rdemmix<-function(nvect,p,g,distr,mu,sigma,dof=NULL,delta=NULL)
{

if(length(c(nvect))!=g) stop("nvect should be a vector")

ndist<-switch(tolower(distr),"mvn"=1,"mvt"=2,"msn"=3,"mst"=4,5)

if(ndist>4) 
stop("the model specified is not available yet")

if(is.null(dof))
dof <- rep(4,g)

if(is.null(delta))
delta <- array(0,c(p,g))

if(length(c(mu)) != (p*g))
stop(paste("mu should be a ",p, 'by', g, "matrix!"))

if(length(c(sigma)) != (p*p*g) )
stop(paste("sigma should be a ",p, 'by', p,'by', g, " array!"))

if(length(c(dof)) != g)
stop(paste("dof should be a ",g, " vector!"))

if(length(c(delta)) != (p*g) )
stop(paste("delta should be a ",p, "by", g, " array!"))

# to fix the "g=1" bug,

mu    = array(mu, c(p,g))
sigma = array(sigma, c(p,p,g))
delta = array(delta, c(p,g))

dat<-array(0,c(10,p))

mvrand<-function(n,p,ndist,mean,cov,nu,del)
{

switch(ndist,
'1' = rdmvn(n,p,mean=mean,cov=cov              ),
'2' = rdmvt(n,p,mean=mean,cov=cov,nu=nu        ),
'3' = rdmsn(n,p,mean=mean,cov=cov,      del=del),
'4' = rdmst(n,p,mean=mean,cov=cov,nu=nu,del=del))
}


if(g>=1)
for(h in 1:g)
{
if(nvect[h]>0)
dat<-rbind(dat,mvrand(nvect[h],p,ndist,mu[,h],sigma[,,h],dof[h],delta[,h]))

}

dat[-(1:10),]
}

