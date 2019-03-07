

#--------------------------------------

# estep procedures

jcamst.estep1<-function(n,p,g,dat,pro,mu,sigma,dof,delta)
{
dat<-as.matrix(dat)
#
tau <- e1 <-e2<-e3<-e4<-e5 <- array(0,c(g*n))
 
obj <- .Fortran('denmst3',PACKAGE="JCM",
as.double(dat),as.integer(n),as.integer(p),as.integer(g),
pro   = as.double(pro),mu    = as.double(mu),sigma = as.double(sigma),
dof   = as.double(dof),delta = as.double(delta),
tau   = as.double(tau),xuu   = as.double(e1),  
elnv = as.double(e2),etw   = as.double(e3),
etw2  = as.double(e4),ewy   = as.double(e5),
loglik= double(1),error = integer(1),as.integer(rep(1,g)))[c(10:17)]

if(obj$error) stop("error")

tau   <- array(obj$tau, c(n,g)) 
elnv <- array(obj$elnv,c(n,g))

sumlnv<- colSums(tau * elnv)

list(error=obj$error,
loglik=     obj$loglik,
tau=  tau,
etw=  array(obj$etw, c(n,g)),
etw2= array(obj$etw2,c(n,g)),
ewy=  array(obj$ewy, c(n,g)),
xuu=  array(obj$xuu, c(n,g)),
sumlnv = sumlnv,sumtau=colSums(tau))
}


jcamst.estep2 <- function(w,etw,etw2,ewy,xuu,y,mu,sigma,dof,delta,theta,thetu,n,p)
{

da <- diag(c(theta),p)
du <- thetu


mmm <- array(0,c(p+1,p+1))
mmm[1:p,1:p]    <- da
mmm[p+1,p+1]    <- du

etee <- array(0,c(p,p))
# 1. 
z1<-diag(mu,p)
z2<-matrix(1,ncol=1,nrow=p)


inv   <-  inverse(sigma,p)

# 2. p by p+1,  w, a, u

xxx <- cbind(z1,z2)

# 3. p+1 by n

# aa + bb * w_i = E(c(a_i,u_i)| y_i, w_i, t_i)

tmp <- t(xxx %*% mmm) %*% inv
aa <- c(rep(1,p),0) + tmp %*% ( t(y) - mu)
bb <- - c(tmp %*% delta) 
cc <- bb %o% (etw*w)

ww <- bb %o% (ewy)

eee <- xxx%*%(aa + ww)+ delta %o% ewy

#eee <- xxx%*%(aa + ww)+ delta %o% ewy

# 4. cov(c(a,u)|y,tau) = E(cov(c(a,u)|w,y,tau)|y,tau)
cc3 <- mmm - tmp %*% (xxx %*% mmm)

# E((y_i - z1 *a - z2 * u) delta^T z_i * w_i | Y_i)
etyz <- ( t(y) - xxx %*% aa ) %*% ( (w*etw) %o% delta) - xxx %*% bb %*% t(delta)* sum(etw2*w)

# E ( E(c(a,u)|y,z) * E(c(a,u)|y,z)|y, z>0)

tmp <- aa %*% ((etw*w) %o% bb)

etauau <- aa %*% (t(aa)*w*xuu) + tmp + t(tmp) + bb %*% t(bb) * sum(etw2*w)

# ystar = y - xxx %*% (a,u)
# E(ystar * t(ystar))

tmp <- (t(y*w*xuu) %*% t(aa) + t(y) %*% t(cc)) %*% t(xxx)
etystar <- t(y) %*% (y*w*xuu) - t(tmp) -tmp + xxx %*% etauau %*% t(xxx)


tmp  <- xxx %*% cc3 %*% t(xxx) *sum(w)
etee <- array(tmp + (etystar - etyz -t(etyz) + delta %*% t(delta) * sum(etw2*w)),c(p,p))

# 8. etee = E(e_i e_i^T | y_i)

#---------------------------------------------------
# eta

# t(t(aa)*w) + cc = E(t_i * c(a_i,u_i| y_i,t_i)

eta <-  matrix((t(t(aa)*w*xuu) + cc )[1:p,],nrow=p)

# dd3 = E( E(t_i*c(a_i,u_i)|y_i,w_i,t_i)  E(c(a_i, u_i,) |w_i,t_i y_i))^T | y_i)

dd3 <- (aa %*% (t(aa)*xuu*w) + aa %*% t(cc) + cc %*% t(aa) + bb %*% t(bb) * sum(etw2*w) )

#------------------------

# 5. etaa 

etaa <- array(cc3[1:p,1:p] * sum(w) + dd3[1:p,1:p],c(p,p))  

# 7. etau

etau <- cc3[1:p,p+1] * sum(w) + dd3[1:p,p+1]

# 8. etuu

etuu <- cc3[p+1,p+1] * sum(w) + dd3[p+1,p+1]


#------------------------


# 6. etaw

if(p==1)
etaw <- sum((aa[1:p,])*etw*w) + bb[1:p] * sum(etw2*w)
else
etaw <- colSums(t(aa[1:p,])*etw*w) + bb[1:p] * sum(etw2*w)


etuwa <- c(colSums( y*etw*w) - xxx %*%( colSums(t(aa)*etw*w) + bb * sum(etw2*w)))

sumtww <- sum(etw2*w)



# 5. etaa2 =E((a_i -1)(a_i -1)^T | y)

tmp <- rowSums(eta)  %o% (rep(1,p)) 

etaa2 <- (etaa - tmp - t(tmp) + matrix(1,ncol=p,nrow=p) * sum(xuu*w) )

etay <- array(0,c(p,p))
for(j in 1:p)
etay[j,] <- colSums(y*eta[j,])

list(etay=etay,etaw=etaw,etuwa=etuwa,etaa=etaa,etau=etau,
etuu=etuu,etee=etee,etaa2=etaa2,sumtww=sumtww,eee=eee)
}


estep.jcamst <- function(datfiles,p,g,pro,mu,sigma,dof,delta,theta,thetu,header,ftype)
{
if((nm <- length(datfiles))<1)
stop("there is no data files in this folder!")

etay <- etaa <- etaa2 <- etee <- array(0,c(p,p,g))

etaw <- etau <- etuwa <- array(0,c(p,g))

sumtww <- sumtau <- sumlnv <- rep(0,g)

etuu <- rep(0,g)

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

#dat  <- as.matrix(read.table(datfiles[j],header=header))
n    <- nrow(dat)

# call library function
eobj <- jcamst.estep1(n,p,g,dat,pro,mu,sigma,dof,delta)

# do each component
for(h in 1:g) 
{

eob2 <- jcamst.estep2(eobj$tau[,h],eobj$etw[,h],eobj$etw2[,h],eobj$ewy[,h],eobj$xuu[,h],dat,
               mu[,h],sigma[,,h],dof[h],delta[,h],theta[,h],thetu[h],n,p)

# sum up 

etay[,,h] <- etay[,,h] + eob2$etay

etaw[,h]  <- etaw[,h]  + eob2$etaw

etuwa[,h] <- etuwa[,h] + eob2$etuwa

etaa[,,h] <- etaa[,,h] + eob2$etaa

etau[,h]  <- etau[,h]  + eob2$etau

etuu[h]   <- etuu[h]   + eob2$etuu

etee[,,h] <- etee[,,h] + eob2$etee

etaa2[,,h]<- etaa2[,,h]+ eob2$etaa2

sumtww[h] <- sumtww[h] + eob2$sumtww

} #end of h loop

#
sumtau <- sumtau + colSums(eobj$tau)

sumlnv<- sumlnv + eobj$sumlnv

loglik <- loglik +  eobj$loglik

nnn <- nnn + n

} #end of j loop

list(loglik=loglik,nnn=nnn,sumtau=sumtau,sumtww=sumtww,sumlnv=sumlnv,
etay=etay,etaw=etaw,etuwa=etuwa,etaa=etaa,etau=etau,
etuu=etuu,etee=etee,etaa2=etaa2)

}


#-------------------------------------
# M step


mstep.jcamst<-function(x,m,g,sumtau,sumtww,sumlnv,sigma,etay,etaw,etuwa,etaa,etau,etuu,etee,etaa2,dofon)
{

alfa <- delta <- theta <- array(0,c(m,g))
thetu<- rep(0,g)

st <- sl <- 0

dof <- rep(3,g)

for( h in 1:g) {

sinv <- inverse(array(sigma[,,h],c(m,m)),m)

if(sumtau[h] >=3)
delta[,h]<- c(etuwa[,h])/sumtww[h]  

# 1. estimate alfa

tmp1 <- rep(0,m)
for(j in 1:m)
tmp1[j] <- sinv[j,]%*% etay[j,,h]

tmp2 <- diag(c(etaw[,h]),m) %*% sinv %*% delta[,h]

tmp3 <- diag(c(etau[,h]),m) %*% sinv %*% rep(1,m)

tmp <- inverse( t(x)%*%(  etaa[,,h] * sinv  )%*%x ,m)

alfa[,h]<-c(tmp %*% ( t(x)%*%(tmp1-tmp2-tmp3)))

# 2. estimate sigma,delta,theta,thetu

if(sumtau[h] >=3)
{
sigma[,,h]<- etee[,,h]/sumtau[h]  

if(m>1) 
theta[,h]<- diag(etaa2[,,h])/sumtau[h]
else
theta[1,h]<- etaa2[1,1,h]/sumtau[h]


thetu[h]<- etuu[h]/sumtau[h]
}

st <- st + sumtau[h]
sl <- sl + sumlnv[h]

if(!dofon)
dof[h]    <- mvt.dof(sumtau[h],sumlnv[h])

} # end h loop

if(dofon)
{
dv <- mvt.dof(st,sl)
for(h in 1:g) dof[h] <- dv
}

list(alfa=alfa,sigma=sigma,dof=dof,delta=delta,theta=theta,thetu=thetu)
}

#---------------------------------------


jcamst.predict <-function(dat,g,pro,mu,sigma,dof,delta,theta,thetu)
{

# input data
dat  <- as.matrix(dat)
n    <- nrow(dat)
p    <- ncol(dat)
mdat<- array(0,c(n,p))
colnames(mdat) <- dimnames(dat)[[2]]


omiga <- array(0,c(p,p,g))
onemat <- array(1,c(p,p))  
#z2 %*% t(z2)

for(h in 1:g) {
z1 <- diag(mu[,h],p)
omiga[,,h] = (sigma[,,h] + array((z1*theta[,h]) %*%t(z1),c(p,p))+ thetu[h] * onemat)
}


# call library function
eobj <- jcamst.estep1(n,p,g,dat,pro,mu,omiga,dof,delta)

clust <- tau2clust(eobj$tau)


# do each component
for(h in 1:g) 
{

eob2 <- jcamst.estep2(eobj$tau[,h],eobj$etw[,h],eobj$etw2[,h],eobj$ewy[,h],eobj$xuu[,h],dat,
               mu[,h],omiga[,,h],dof[h],delta[,h],theta[,h],thetu[h],n,p)

mdat<- (mdat + t(eob2$eee) * ifelse(clust==h,1,0))


} #end of h loop

list(eee=mdat,clust=clust,tau=eobj$tau)

}

jcamst <-function(datfiles,p,g,initobj,ncov=3,dofon=0,debug=1,
itmax=1000,epsilon=1e-5,header=TRUE,ftype="csv" )
{

x  <- diag(p)

alfa     <- initobj$mu
sigma    <- initobj$sigma

theta    <- initobj$theta
thetu    <- initobj$thetu
delta    <- initobj$delta 
dof      <- initobj$dof
pro      <- initobj$pro


onemat <- array(1,c(p,p))  #z2 %*% t(z2)


mu <- array(0,c(p,g))
omiga <- array(0,c(p,p,g))


# controls
flag<-0;lk<-rep(0,itmax)
#
for(it in 1:itmax)
{
# e-step

for(h in 1:g) {
mu[,h] <-c( x%*%alfa[,h])
z1 <- diag(mu[,h],p)

omiga[,,h] = sigma[,,h] + array((z1*theta[,h]) %*%t(z1),c(p,p))+ thetu[h] * onemat
}

ee <- estep.jcamst(datfiles,p,g,pro,mu,omiga,dof,delta,theta,thetu,header,ftype)
pro   <- ee$sumtau/ee$nnn


lk[it]<-ee$loglik
cat('\n', it,lk[it])

mobj <- mstep.jcamst(x,p,g,ee$sumtau,ee$sumtww,ee$sumlnv,
sigma,ee$etay,ee$etaw,ee$etuwa,ee$etaa,ee$etau,
ee$etuu,ee$etee,ee$etaa2,dofon)


# m-step

delta <- mobj$delta
alfa  <- mobj$alfa
sigma <- mobj$sigma

theta <- mobj$theta
thetu <- mobj$thetu

dof   <- mobj$dof

if(ncov!=3)
sigma <-  getcov(sigma,ee$sumtau,ee$nnn,p,g,ncov)


loglik <- lk[it]

if(it<11) next

if(abs(lk[it]-lk[it-10])<epsilon*abs(lk[it-10]) ) {flag<-1;break}
}


n <- ee$nnn

# bic

nk <- switch(ncov,'1'=p*(p+1)/2,'2'=p,'3'= p*(p+1)/2*g,'4'=p*g)
nk <- nk + (g-1) + 3*p*g + 2*g
bic <- -2*loglik + log(n)*nk


list(error=1-flag,loglik=loglik,bic=bic,pro=pro,alfa=alfa,mu=mu, 
sigma=sigma,dof=dof,delta=delta,theta=theta,thetu=thetu,distr="mst",class="jcm")
}

#-----------------------------------------------

