
# E step
msmvn.estep1<-function(y,n,p,g,pro,mu,sigma)
{

obj <- .Fortran('estepmvn',PACKAGE="JCM",
as.double(y),as.integer(n),as.integer(p),as.integer(g),
as.double(pro),as.double(mu),as.double(sigma),
tau=double(n*g),sumtau = double(g),
loglik= double(1),error = integer(1))[8:11]

if(obj$error) stop("error")
tau   <- array(obj$tau, c(n,g)) 
list(tau=tau,sumtau=obj$sumtau,loglik=obj$loglik)
}



msmvn.estep2 <- function(y,n,p,g,tau,mu)
{

obj <- .Fortran('scaestepmvn',PACKAGE="JCM",
as.double(y),as.integer(n),as.integer(p),as.integer(g),
as.double(tau),as.double(mu),
ewy = double(p*g),ewyy = double(p*p*g))[7:8]

list(ewy=array(obj$ewy,c(p,g)),ewyy=array(obj$ewyy,c(p,p,g)))

}

# do multiple samples
estep.msmvn <- function(datfiles,p,g,pro,mu,sigma,header,ftype)
{
if((nm <- length(datfiles))<1)
stop("there is no data files in this folder!")

ewy <- array(0,c(p,g))

sumtau <- rep(0,g)

ewyy <- array(0,c(p,p,g))

loglik <- 0

nnn <- 0


for(j in 1:nm)
{
# input data
if(tolower(ftype)=="txt")
  dat  <- as.matrix(read.table(datfiles[j],header=header))[,1:p]

if(tolower(ftype)=="csv")
  dat  <- as.matrix(read.csv(datfiles[j],header=header))[,1:p]


n <- nrow(dat)

#esteps

eobj  <- msmvn.estep1(dat,n,p,g,pro,mu,sigma)

eob2  <- msmvn.estep2(dat,n,p,g,eobj$tau,mu)

# do each component
for(h in 1:g) 
{
# sum up 
ewyy[,,h] <- ewyy[,,h] + eob2$ewyy[,,h]
ewy[,h]   <- ewy[,h]   + eob2$ewy[,h]
} #end of h loop

#
sumtau <- sumtau +  eobj$sumtau

loglik <- loglik +  eobj$loglik

nnn <- nnn + n

} #end of j loop

list(loglik=loglik,nnn=nnn,sumtau=sumtau,ewy=ewy,ewyy=ewyy)

}
#end fun


# M step
mstep.msmvn<-function(p,g,sumtau,ewy,ewyy)
{

mu    <- array(0,c(p,g)  )
sigma <- array(0,c(p,p,g))

for( h in 1:g) {

# estimate mu, sigma

if(sumtau[h] >=3)
{
sigma[,,h] <- ewyy[,,h]/sumtau[h]   
mu[,h]     <- ewy[,  h]/sumtau[h]
}

} # end h loop

list(mu=mu,sigma=sigma)
}



#  Fit the Mixture models
msmvn<-function(datfiles,p,g,initobj,ncov=3,debug=1,
itmax=1000,epsilon=1e-5,header=TRUE,ftype="txt")
{
#initials
mu       <- initobj$mu
sigma    <- initobj$sigma
pro      <- initobj$pro

# controls
flag<-0;lk<-rep(0,itmax)
#
for(it in 1:itmax)
{
# e-step
ee    <- estep.msmvn(datfiles,p,g,pro,mu,sigma,header,ftype)

pro   <- ee$sumtau/ee$nnn

lk[it]<- ee$loglik

if(debug)
cat('\n', it,lk[it])
mobj <- mstep.msmvn(p,g,ee$sumtau,ee$ewy,ee$ewyy)
# m-step
mu    <- mobj$mu
sigma <- mobj$sigma
if(ncov!=3)
sigma <-  getcov(sigma,ee$sumtau,ee$nnn,p,g,ncov)
loglik <- lk[it]
if(it<11) next
if(abs(lk[it]-lk[it-10])<epsilon*abs(lk[it-10]) ) {flag<-1;break}
}

list(error=1-flag,loglik=loglik,nnn=ee$nnn,
pro=pro,mu=mu,sigma=sigma,dof=rep(0,g),delta=array(0,c(p,g)),distr="mvn")
}
#end of program


