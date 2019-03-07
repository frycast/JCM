

#--------------------------------------

# estep procedures



jcamvt.estep1<-function(n,p,g,dat,pro,mu,sigma,dof)
{
dat<- as.matrix(dat)

tau <- xuu  <- array(0,c(n,g))
 
obj <- .C('estepmvt_',PACKAGE="JCM",
as.double(dat),
as.integer(n),as.integer(p),as.integer(g),
pro   = as.double(pro),
mu    = as.double(mu),
sigma = as.double(sigma),
dof   = as.double(dof),
tau   = as.double(tau),
xuu   = as.double(xuu),
sumtau= double(g),
sumxuu= double(g),
sumlnv=double(g),
loglik= double(1),error = integer(1) )[c(9:15)]

if(obj$error) stop("error")

list(error=obj$error,
loglik=     obj$loglik,
tau=  array(obj$tau, c(n,g)),
xuu=  array(obj$xuu, c(n,g)),
sumtau=c(obj$sumtau),sumlnv=c(obj$sumlnv))
}

# for h component

jcamvt.estep2<-function(n,p,y,w,xuu,mu,omiga,dof,theta,thetu)
{
#theta is a vector

acov<-diag(theta)

mmm<-matrix(0,ncol=p+1,nrow=p+1)
mmm[1:p,1:p]    <- acov
mmm[p+1,p+1]    <- thetu

#---------------------------

z1 <- diag(c(mu))
z2 <- matrix(1,ncol=1,nrow=p)
#inv<-solve(omiga)
inv<-inverse(omiga,p)

#---------------------------

# E(a|y),cov(a|y)
eea<-rep(1,p)+acov%*%t(z1)%*%inv%*%(t(y)-mu)
eca<-(acov-acov%*%t(z1)%*%inv%*%z1%*%acov)  


# E(u|y), cov(u|y)
eeu<-thetu * t(z2)%*%inv%*%(t(y)-mu)
ecu<-(thetu-thetu^2 * t(z2)%*%inv%*%z2)

eee1<-t(y)-z1%*%eea-z2%*%eeu

ece<-(cbind(z1,z2)%*%(mmm-rbind(acov%*%t(z1),thetu * t(z2))%*%inv%*%cbind(z1%*%t(acov),z2 * thetu ))%*%rbind(t(z1),t(z2))) 


# E(A,U|Y)

ecau <- (-t(acov)%*%t(z1)%*%inv%*%z2 * thetu )


#--------------------------

#do alpha

etay <- array(0,c(p,p))
for(j in 1:p)
etay[j,]  <- colSums( y *      (eea[j,]*w*xuu)) 

etau     <- (ecau*sum(w)+colSums(t(eea)* c(eeu)*w*xuu ))

etaa1    <- (eca*sum(w) + eea%*%(t(eea)*w*xuu ))

#--------------------------

# do sigma, theta, thetau

etee     <- (ece*sum(w) +  eee1   %*%(t(eee1)  *w*xuu))

etaa2    <- (eca*sum(w) + (eea-rep(1,p))%*%(t(eea-rep(1,p))*w*xuu))

etuu     <- (ecu*sum(w) +  eeu   %*%(t(eeu)  *w*xuu))

#--------------------------
eee <- z1%*%eea + z2%*%eeu

list(etay=etay,etau=etau,etaa1=etaa1,etee=etee,etaa2=etaa2,etuu=etuu,eee=eee)

}

#-------------------------------

estep.jcamvt <- function(datfiles,p,g,pro,mu,sigma,dof,theta,thetu,header,ftype)
{
if((nm <- length(datfiles))<1)
stop("there is no data files in this folder!")

etay <- etaa1 <- etaa2 <- etee <- array(0,c(p,p,g))

etau <- array(0,c(p,g))

sumtau <-  sumlnv <- rep(0,g)

etuu   <- rep(0,g)

loglik <- 0

nnn <- 0


for(j in 1:nm)
{

# input data

if(tolower(ftype)=="csv")
dat <- read.csv(datfiles[j],header=header)[,1:p]

if(tolower(ftype)=="txt")
dat <- read.table(datfiles[j],header=header)[,1:p]

dat  <- as.matrix(dat)
n    <- nrow(dat)

# call library function
eobj <- jcamvt.estep1(n,p,g,dat,pro,mu,sigma,dof)

# do each component
for(h in 1:g) 
{

eob2 <- jcamvt.estep2(n,p,dat,eobj$tau[,h],eobj$xuu[,h],
               mu[,h],sigma[,,h],dof[h],theta[,h],thetu[h])

# sum up 

etay[,,h]  <- etay[,,h] + eob2$etay

etaa1[,,h] <- etaa1[,,h] + eob2$etaa1

etau[,h]  <- etau[,h]  + eob2$etau

etuu[h]   <- etuu[h]   + eob2$etuu

etee[,,h] <- etee[,,h] + eob2$etee

etaa2[,,h]<- etaa2[,,h]+ eob2$etaa2


} #end of h loop

#
sumtau <- sumtau +  eobj$sumtau

loglik <- loglik +  eobj$loglik

sumlnv<- sumlnv + eobj$sumlnv


nnn <- nnn + n

} #end of j loop

list(loglik=loglik,nnn=nnn,sumtau=sumtau,sumlnv=sumlnv,
etay=etay,etau=etau,etaa1=etaa1,etee=etee,etaa2=etaa2,etuu=etuu)

}
#end fun


#-------------------------------------
# M step


mstep.jcamvt <-function(p,g,x,sigma,sumtau,sumlnv,etay,etau,etaa1,etee,etaa2,etuu,model)
{

alfa <- theta <- array(0,c(p,g))
thetu<- rep(0,g)
z2 <- matrix(1,ncol=1,nrow=p)

sl <- st <- 0

dof <- rep(3,g)

for( h in 1:g) {

sinv <- inverse(array(sigma[,,h],c(p,p)),p)


# 1. estimate alfa

tmp1 <- rep(0,p)
for(j in 1:p)
tmp1[j] <- sinv[j,]%*% etay[j,,h]

tmp2 <- (sinv%*%z2)*etau[,h]

tmp <- inverse( t(x)%*%(  etaa1[,,h] * sinv  )%*%x ,p)

alfa[,h]<-c(tmp %*% ( t(x)%*%(tmp1-tmp2)))

# 2. estimate sigma,delta,theta,thetu

if(sumtau[h] >=3)
{
sigma[,,h]<- etee[,,h]/sumtau[h]  

if(p>1) 
theta[,h]<- diag(etaa2[,,h])/sumtau[h]
else
theta[1,h]<- etaa2[1,1,h]/sumtau[h]


thetu[h]<- etuu[h]/sumtau[h]
}

if(!model)
dof[h]    <- mvt.dof(sumtau[h],sumlnv[h])

st = st + sumtau[h]
sl = sl + sumlnv[h]

} # end h loop

if(model)
{
dv    <- mvt.dof(st,sl)
for(h in 1:g)
dof[h]<- dv
}

list(alfa=alfa,sigma=sigma,dof=dof,theta=theta,thetu=thetu)
}

#-----------------------------------------------

jcamvt.predict <-function(dat,g,pro,mu,sigma,dof,theta,thetu)
{

# input data
dat  <- cbind(dat)
dat<- as.matrix(dat)


n    <- nrow(dat)
p    <- ncol(dat)
mdat<- array(0,c(n,p))
colnames(mdat) <- dimnames(dat)[[2]]


omiga <- array(0,c(p,p,g))
onemat <- array(1,c(p,p))  #z2 %*% t(z2)

for(h in 1:g) {
z1 <- diag(mu[,h],p)
omiga[,,h] = (sigma[,,h] + array((z1*theta[,h]) %*%t(z1),c(p,p))+ thetu[h] * onemat)
}


eobj <- jcamvt.estep1(n,p,g,dat,pro,mu,omiga,dof)
clust <- tau2clust(eobj$tau)

# do each component
for(h in 1:g) 
{

eob2 <- jcamvt.estep2(n,p,dat,eobj$tau[,h],eobj$xuu[,h],
               mu[,h],omiga[,,h],dof[h],theta[,h],thetu[h])

# prediction
mdat<- (mdat + t(eob2$eee) * ifelse(clust==h,1,0))


} #end of h loop

list(eee=mdat,clust=clust,tau=eobj$tau)
}


jcamvt <-function(datfiles,p,g,initobj,ncov=3,dofon=0,debug=1,
itmax=1000,epsilon=1e-5,header=TRUE,ftype="csv" )
{

x  <- diag(p)

#initials

pro      <- initobj$pro
alfa     <- initobj$mu
sigma    <- initobj$sigma
dof      <- initobj$dof

theta    <- initobj$theta
thetu    <- initobj$thetu


onemat <- array(1,c(p,p))  #z2 %*% t(z2)


mu    <- array(0,c(p,g))
omiga <- array(0,c(p,p,g))


# controls
flag<-0;lk<-rep(0,itmax)
#
for(it in 1:itmax)
{

# e-step
#-------------------------

# 1. get mu and sigma

for(h in 1:g) {
mu[,h] <-c( x%*%alfa[,h])
z1 <- diag(mu[,h],p)
omiga[,,h] = (sigma[,,h] + array((z1*theta[,h]) %*%t(z1),c(p,p))+ thetu[h] * onemat)
}

# 2. call estep

ee   <- estep.jcamvt(datfiles,p,g,pro,mu,omiga,dof,theta,thetu,header,ftype)
pro  <- ee$sumtau/ee$nnn


loglik  <-ee$loglik
cat('\n', it,loglik)

# 3. call mstep

mobj <- mstep.jcamvt(p,g,x,sigma,ee$sumtau,ee$sumlnv,
ee$etay,ee$etau,ee$etaa1,ee$etee,ee$etaa2,ee$etuu,dofon)


# 4. update the parameters

alfa  <- mobj$alfa
sigma <- mobj$sigma

theta <- mobj$theta
thetu <- mobj$thetu

dof   <- mobj$dof

if(ncov!=3)
sigma <-  getcov(sigma,ee$sumtau,ee$nnn,p,g,ncov)


lk[it] <- loglik 

if(it<11) next

if(abs(lk[it]-lk[it-10])<epsilon*abs(lk[it-10]) ) {flag<-1;break}
}


n <- ee$nnn

# bic

nk <- switch(ncov,'1'=p*(p+1)/2,'2'=p,'3'= p*(p+1)/2*g,'4'=p*g)
nk <- nk + (g-1) + 2*p*g + 2*g
bic <- -2*loglik + log(n)*nk


list(error=1-flag,loglik=loglik,bic=bic,pro=pro,alfa=alfa,mu=mu, 
sigma=sigma,dof=dof,delta=array(0,c(p,g)),theta=theta,thetu=thetu,
distr="mvt",class="jcm")
}

#--------------------------

