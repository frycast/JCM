#  Pooling Functions 
msmvt.estep1<-function(dat,n,p,g,pro,mu,sigma,dof)
{

obj <- .Fortran('estepmvt',PACKAGE="JCM",
as.double(dat),as.integer(n),as.integer(p),as.integer(g),
as.double(pro),as.double(mu),as.double(sigma),as.double(dof),
tau= double(n*g),xuu= double(n*g),
sumtau = double(g),sumxuu = double(g),sumxuuln = double(g),
loglik= double(1),error = integer(1))[-(1:8)]

if(obj$error) stop("error")

list(tau=array(obj$tau, c(n,g)) ,xuu=array(obj$xuu,c(n,g)),
sumtau=obj$sumtau,sumxuu=obj$sumxuu,
sumxuuln=obj$sumxuuln,loglik=obj$loglik)
}

msmvt.estep2 <- function(y,n,p,g,tau,xuu,mu)
{

obj <- .Fortran('scaestepmvt',PACKAGE="JCM",
as.double(y),as.integer(n),as.integer(p),as.integer(g),
as.double(tau),as.double(xuu),as.double(mu),
ewy = double(p*g),ewyy = double(p*p*g))[-(1:7)]

list(ewy=array(obj$ewy,c(p,g)),ewyy=array(obj$ewyy,c(p,p,g)))
}


estep.msmvt <- function(datfiles,p,g,pro,mu,sigma,dof,header,ftype)
{
if((nm <- length(datfiles))<1)
stop("there is no data files in this folder!")
ewyy <- array(0,c(p,p,g))
ewy <- array(0,c(p,g))
sumtau <- sumxuu <- sumxuuln <- rep(0,g)
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

eobj <- msmvt.estep1(dat,n,p,g,pro,mu,sigma,dof)

eob2 <- msmvt.estep2(dat,n,p,g,eobj$tau,eobj$xuu,mu)

# do each component
for(h in 1:g) 
{
# sum up 
ewy[,h]   <- ewy[,h]   + eob2$ewy[,h]
ewyy[,,h] <- ewyy[,,h] + eob2$ewyy[,,h]

} #end of h loop

#tau and eee and clust are written to data files
#
sumtau <- sumtau + eobj$sumtau
loglik <- loglik + eobj$loglik
sumxuu <- sumxuu + eobj$sumxuu
sumxuuln <- sumxuuln + eobj$sumxuuln
nnn    <- nnn + n
} #end of j loop

list(loglik=loglik,nnn=nnn,sumtau=sumtau,sumxuuln=sumxuuln,sumxuu=sumxuu,
ewy=ewy,ewyy=ewyy)

}
#end fun


# M step
mstep.msmvt<-function(p,g,sumtau,sumxuu,sumxuuln,ewy,ewyy,model=0)
{

st<- sl <- 0
dof    <- rep(50,g)
mu     <- array(0,c(p,g))
sigma  <- array(0,c(p,p,g))
for( h in 1:g) {
if(sumtau[h] >=3)
{
mu[,h]<- ewy[,h]/sumxuu[h]
sigma[,,h]<- ewyy[,,h]/sumtau[h]  
}

if(!model)
dof[h]    <- mvt.dof(sumtau[h],sumxuuln[h])
st = st + sumtau[h]
sl = sl + sumxuuln[h]
} # end h loop
if(model)
{
dv    <- mvt.dof(st,sl)
for(h in 1:g)
dof[h]<- dv
}
list(mu=mu,sigma=sigma,dof=dof)
}

#-------------------------------------

#  Fit the Mixture models
msmvt<-function(datfiles,p,g,initobj,ncov=3,common=0,
debug=1,itmax=1000,epsilon=1e-5,header=TRUE,ftype="txt" )
{

#initials

pro      <- initobj$pro
mu       <- initobj$mu
sigma    <- initobj$sigma
dof      <- initobj$dof

# controls
flag<-0;lk<-rep(0,itmax)
#
for(it in 1:itmax)
{
# e-step
#-------------------------
# call estep
ee   <- estep.msmvt(datfiles,p,g,pro,mu,sigma,dof,header,ftype)

pro  <- ee$sumtau/ee$nnn

loglik  <-ee$loglik
cat('\n', it,loglik)
#  call mstep
mobj <- mstep.msmvt(p,g,ee$sumtau,ee$sumxuu,ee$sumxuuln,ee$ewy,ee$ewyy,common)
# update the parameters
mu    <- mobj$mu
sigma <- mobj$sigma
dof   <- mobj$dof
if(ncov!=3)
sigma <-  getcov(sigma,ee$sumtau,ee$nnn,p,g,ncov)
lk[it] <- loglik 
if(it<11) next
if(abs(lk[it]-lk[it-10])<epsilon*abs(lk[it-10]) ) {flag<-1;break}
}
n <- ee$nnn
# bic
list(error=1-flag,loglik=loglik,nnn=n,pro=pro,mu=mu, 
sigma=sigma,dof=dof,delta=array(0,c(p,g)),distr="mvt")
}

# end of code
