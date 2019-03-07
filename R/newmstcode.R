#  Pooling Functions 

msmst.estep1<-function(dat,n,p,g,pro,mu,sigma,dof,delta)
{
#
obj <- .Fortran('estepmst',PACKAGE="JCM",
as.double(dat),as.integer(n),as.integer(p),as.integer(g),
as.double(pro),as.double(mu),as.double(sigma),as.double(dof),as.double(delta),
tau= double(n*g),ev= double(n*g),elnv= double(n*g),
ez1v= double(n*g),ez2v= double(n*g),
sumtau = double(g),sumvt = double(g),sumzt = double(g),sumlnv = double(g),
loglik= double(1),error = integer(1),as.integer(rep(1,g)))[-(1:9)]

if(obj$error) stop("error")

list(tau=array(obj$tau, c(n,g)) ,ev=array(obj$ev,c(n,g)),
elnv=array(obj$elnv,c(n,g)),ez1v=array(obj$ez1v,c(n,g)),
ez2v=array(obj$ez2v,c(n,g)),
sumtau=obj$sumtau,sumvt=obj$sumvt,sumzt=obj$sumzt,sumlnv=obj$sumlnv,
loglik=obj$loglik)
}



msmst.estep2 <- function(y,n,p,g,tau,ev,ez1v,ez2v,mu,delta)
{

obj <- .Fortran('scaestepmst',PACKAGE="JCM",
as.double(y),as.integer(n),as.integer(p),as.integer(g),
as.double(tau),as.double(ev),
as.double(ez1v),as.double(ez2v),
as.double(mu),as.double(delta),
ewy = double(p*g),ewz = double(p*g),ewyy = double(p*p*g))[-(1:10)]

list(ewy=array(obj$ewy,c(p,g)),ewz=array(obj$ewz,c(p,g)),ewyy=array(obj$ewyy,c(p,p,g)))
}

estep.msmst <- function(datfiles,p,g,pro,mu,sigma,dof,delta,header,ftype)
{
if((nm <- length(datfiles))<1)
stop("there is no data files in this folder!")

ewyy <- array(0,c(p,p,g))
ewy <- ewz <- array(0,c(p,g))
 
sumtau <- sumvt <-sumzt <-sumlnv <- rep(0,g)

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

eobj <- msmst.estep1(dat,n,p,g,pro,mu,sigma,dof,delta)
clust <- tau2clust(eobj$tau)

mdat<- array(0,c(n,p))

eob2 <- msmst.estep2(dat,n,p,g,eobj$tau,eobj$ev,eobj$ez1v,eobj$ez2v,mu,delta)

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
sumvt  <- sumvt  +  eobj$sumvt
sumzt  <- sumzt  +  eobj$sumzt
sumlnv <- sumlnv +  eobj$sumlnv
loglik <- loglik +  eobj$loglik
nnn <- nnn + n
} #end of j loop

list(loglik=loglik,nnn=nnn,sumtau=sumtau,sumvt=sumvt,
sumzt=sumzt,sumlnv=sumlnv,
ewy=ewy,ewz=ewz,ewyy=ewyy)

}
#end fun


# M step
mstep.msmst<-function(p,g,sumtau,sumvt,sumzt,sumlnv,ewy,ewz,ewyy)
{

mu     <- delta    <-array(0,c(p,g))
sigma  <- array(0,c(p,p,g))
dof<-rep(0,g)

for( h in 1:g) {

if(sumtau[h] >=3)
{
mu[,h]<- ewy[,h]/sumvt[h]
delta[,h]<- ewz[,h]/sumzt[h]
sigma[,,h]<- ewyy[,,h]/sumtau[h]  
}

dof[h]<- mvt.dof(sumtau[h],sumlnv[h])


}


list(mu=mu,sigma=sigma,dof=dof,delta=delta)
}

#-------------------------------------

#  Fit the Mixture models
msmst<-function(datfiles,p,g,initobj,ncov=3,
debug=1,itmax=1000,epsilon=1e-5,header=TRUE,ftype="txt" )
{

#initials

pro      <- initobj$pro
mu       <- initobj$mu
sigma    <- initobj$sigma
delta    <- initobj$delta
dof      <- initobj$dof

# controls
flag<-0;lk<-rep(0,itmax)
#
for(it in 1:itmax)
{
# e-step
#-------------------------
# call estep
ee   <- estep.msmst(datfiles,p,g,pro,mu,sigma,dof,delta,header,ftype)

pro  <- ee$sumtau/ee$nnn

loglik  <-ee$loglik
cat('\n', it,loglik)

#  call mstep
mobj <- mstep.msmst(p,g,ee$sumtau,ee$sumvt,ee$sumzt,ee$sumlnv,ee$ewy,ee$ewz,ee$ewyy)

# update the parameters
mu    <- mobj$mu
sigma <- mobj$sigma
delta <- mobj$delta
dof   <- mobj$dof

if(ncov!=3)
sigma <-  getcov(sigma,ee$sumtau,ee$nnn,p,g,ncov)
lk[it] <- loglik 
if(it<11) next
if(abs(lk[it]-lk[it-10])<epsilon*abs(lk[it-10]) ) {flag<-1;break}
}
n <- ee$nnn

list(error=1-flag,loglik=loglik,nnn=n,pro=pro,mu=mu, 
sigma=sigma,delta=delta,dof=dof,distr="mst")
}

# end of code
