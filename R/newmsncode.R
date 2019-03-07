#  Pooling Functions 
msmsn.estep1<-function(y,n,p,g,pro,mu,sigma,delta)
{

obj <- .Fortran('estepmsn',PACKAGE="JCM",
as.double(y),as.integer(n),as.integer(p),as.integer(g),
as.double(pro),as.double(mu),as.double(sigma),as.double(delta),
tau= double(n*g),ev= double(n*g),vv= double(n*g),
sumtau = double(g),sumev = double(g),
loglik= double(1),error = integer(1))[9:15]

if(obj$error) stop("error")

list(tau=array(obj$tau,c(n,g)),ev=array(obj$ev,c(n,g)),
vv=array(obj$vv,c(n,g)),sumtau=obj$sumtau,sumev=obj$sumev,
loglik=obj$loglik)
}

msmsn.estep2 <- function(y,n,p,g,tau,ev,vv,mu,delta)
{
obj <- .Fortran('scaestepmsn',PACKAGE="JCM",
as.double(y),as.integer(n),as.integer(p),as.integer(g),
as.double(tau),as.double(ev),as.double(vv),
as.double(mu),as.double(delta),
ewy = double(p*g),ewz = double(p*g),ewyy = double(p*p*g))[10:12]

list(ewy=array(obj$ewy,c(p,g)),ewz=array(obj$ewz,c(p,g)),
ewyy=array(obj$ewyy,c(p,p,g)))

}

estep.msmsn <- function(datfiles,p,g,pro,mu,sigma,delta,header,ftype)
{
if((nm <- length(datfiles))<1)
stop("there is no data files in this folder!")

ewyy <- array(0,c(p,p,g))

ewy <- ewz <- array(0,c(p,g))
 
sumtau <- sumev <- rep(0,g)

loglik <- 0

nnn <- 0

for(j in 1:nm)
{
# input data
if(tolower(ftype)=="txt")
  dat  <- as.matrix(read.table(datfiles[j],header=header))[,1:p]

if(tolower(ftype)=="csv")
  dat  <- as.matrix(read.csv(datfiles[j],header=header))[,1:p]


n    <- nrow(dat)

eobj <- msmsn.estep1(dat,n,p,g,pro,mu,sigma,delta)

eob2 <- msmsn.estep2(dat,n,p,g,eobj$tau,eobj$ev,eobj$vv,mu,delta)

# do each component
for(h in 1:g) 
{
# sum up 
ewy[,h]   <- ewy[,h]   + eob2$ewy[,h]
ewz[,h]   <- ewz[,h]   + eob2$ewz[,h]
ewyy[,,h] <- ewyy[,,h] + eob2$ewyy[,,h]

} #end of h loop

#tau and eee and clust are written to data files
#
sumtau <- sumtau +  eobj$sumtau
sumev  <- sumev  +  eobj$sumev
loglik <- loglik +  eobj$loglik
nnn <- nnn + n
} #end of j loop

list(loglik=loglik,nnn=nnn,sumtau=sumtau,sumev=sumev,
ewy=ewy,ewz=ewz,ewyy=ewyy)

}
#end fun


# M step
mstep.msmsn<-function(p,g,sumtau,sumev,ewy,ewz,ewyy)
{

mu     <- delta    <-array(0,c(p,g))
sigma  <- array(0,c(p,p,g))

for( h in 1:g) {
if(sumtau[h] >=3)
{
mu[,h]    <- ewy[,h]/sumtau[h]
delta[,h] <- ewz[,h]/sumev[h]
sigma[,,h]<- ewyy[,,h]/sumtau[h]  
}

}
list(mu=mu,sigma=sigma,delta=delta)
}

#-------------------------------------

#  Fit the Mixture models
msmsn<-function(datfiles,p,g,initobj,ncov=3,
debug=1,itmax=1000,epsilon=1e-5,header=TRUE,ftype="txt" )
{

#initials

pro      <- initobj$pro
mu       <- initobj$mu
sigma    <- initobj$sigma
delta    <- initobj$delta

# controls
flag<-0;lk<-rep(0,itmax)
#
for(it in 1:itmax)
{
# e-step
#-------------------------
# call estep
ee   <- estep.msmsn(datfiles,p,g,pro,mu,sigma,delta,header,ftype)

pro  <- ee$sumtau/ee$nnn

loglik  <-ee$loglik
cat('\n', it,loglik)
#  call mstep
mobj <- mstep.msmsn(p,g,ee$sumtau,ee$sumev,ee$ewy,ee$ewz,ee$ewyy)

# update the parameters
mu    <- mobj$mu
sigma <- mobj$sigma
delta   <- mobj$delta

if(ncov!=3)
sigma <-  getcov(sigma,ee$sumtau,ee$nnn,p,g,ncov)
lk[it] <- loglik 
if(it<11) next
if(abs(lk[it]-lk[it-10])<epsilon*abs(lk[it-10]) ) {flag<-1;break}
}
n <- ee$nnn

list(error=1-flag,loglik=loglik,nnn=n,pro=pro,mu=mu, 
sigma=sigma,delta=delta,dof=rep(0,g),distr="msn")
}

# end of code
