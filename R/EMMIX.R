#*****************************************************
#
#  EM algorithm for Mixture of Multivariate (Skew) Normal/T Distributioins
#
#  Package: EMMIX-JCM
#
#  Version: 1.0-20
#
#  Code by Kui (Sam) Wang <kwang@maths.uq.edu.au>
#  
#*****************************************************
#  
#  Changes:
#  
#  Version 1.0-19:

#  1) A bug was fixed in function emmix, when initial cluster labels are 
#     not in the range of 1 to g, error occurs. 
#     checks on input of "clust" and "distr" are added. 
#  2) Supervised Clustering Methods are incorporated.
#  3) Discriminant Analysis via Kullback-Leibler distance (KLD) is incorporated.
#  4) Some pairwise plots are added to plot the heatmap of FCS data and 
#     contours of the fitted mixture models.
#  5) The function emmixMOD is removed and replaced by emmixMOD2.

#  Version < 1.0-18
#
#  1) The old versions called "emmixskew" are no longer upated.
#  2) UNIX version "EMMIX 1.3" on our website will not be updated.
#  3) R package EMMIX is a majoy update for above two projects.

#  Version 1.0-18:

#  1) A bug was found by Michael Fahey <Michael.Fahey@csiro.au>
#  errors occur when g=1.
# 
#  2) emmixMOD2() is added for speeding up 
#     the calculations of mode point for each component.


#  
#****************************************************** 

.packageName <- 'JCM'

.First.lib<-function(lib,pkg) { library.dynam('JCM',pkg,lib) }


emmix<-function(dat, g, distr = "mvn",ncov= 3,   
clust   = NULL,init= NULL,  
itmax   = 1000,epsilon = 1e-6,nkmeans = 0, nrandom = 10, 
nhclust = FALSE,debug = TRUE,initloop=20)
{
dat<- as.matrix(dat)

if(!is.null(init)  | !missing(init))
obj<-emmixfit2(dat,g, init, distr,ncov,itmax,epsilon)
else
{
if(is.null(clust) | missing(clust)) {

init<- try(init.mix(dat,g,distr,ncov,nkmeans,nrandom,nhclust,initloop))

if(!is.null(init)) 
obj<-emmixfit2(dat,g, init, distr,ncov,itmax,epsilon)
else
{
warning("not find initial values")
obj<-list()
obj$error=20
}
}
else {
obj<-emmixfit1(dat,g,clust, distr,ncov,itmax,epsilon,initloop)
}
}

error <- obj$error;

ret<-NULL

msg <- switch(tolower(error),
'1' = paste("stopped at (did not converge within) ", itmax, " iterations"),
'2' = paste("density fails at initial steps! "),
'3' = paste("allocation fails at initial steps"),
'12' = paste("density fails at estps! "),
'13' = paste("allocation fails at esteps"),
'20' =paste("not find initials"))


if(error>1)
{
cat('\n-----------------------\n')
warning("error code=",error,'\n',msg,"\n")
cat('\n-----------------------\n\n')
}

# summarize the results

if(error<=1) {

ICL <- getICL(dat,nrow(dat),ncol(dat),g,distr,ncov,obj$pro,obj$mu,obj$sigma,obj$dof,obj$delta,obj$clust)

ret<-obj

ret$ICL<-ICL$ICL

# mode point for each component

ret$modpts<- emmixMOD(ncol(dat),g,distr,ret$mu,ret$sigma,ret$dof,ret$delta)


if(debug) {

msg<-switch(tolower(distr), 
'mvn'=paste(g,"- Component Multivariate Normal Mixture Model"),
'mvt'=paste(g,"- Component Multivariate t      Mixture Model"),
'msn'=paste(g,"- Component Multivariate Skew Normal Mixture Model"),
'mst'=paste(g,"- Component Multivariate Skew-t Mixture Model"))


cat('\n-----------------------\n\n')
cat(msg,"\n")
cat('\n-----------------------\n\n')
 
switch(tolower(distr), 'mvn' = print(obj[1:8]),'mvt' = print(obj[1:9]),
'msn' = print(obj[c(1:8,10)]),'mst' = print(obj[1:10]))
print(ICL)
cat('\n-----------------------\n')
}

}
ret
}

getICL<-function(x,n,p,g,distr,ncov,pro,mu,sigma,dof,delta,clust)
{

x<-as.matrix(x);loglik<-0;nc=0

#if(length(table(clust<- unclass(as.ordered(clust))))!=g)
#stop(paste("labels should be of ", g, "levels"))

nn <- sum(clust>0) #outliers are numbered zero.

lnden<-(as.matrix((ddmix(x,n,p,g, distr, mu, sigma, dof, delta))))

for(h in 1:g) loglik=loglik+sum(ifelse(clust==h,log(pro[h])+lnden[,h],0) )

nc<-switch(tolower(ncov),
    '1' = (g-1)+g*p+p*(1+p)/2,  #common covariance
    '2' = (g-1)+g*p+p,          #common diagonal covariance
    '3' = (g-1)+g*(p+p*(1+p)/2),#general covariance
    '4' = (g-1)+g*(p+p),        #general diagonal covariance
    '5' = (g-1)+g*(p+1)  )      #sigma(h)*I_p

nc <- switch(tolower(distr),"mvn" = nc,"mvt" = nc + g,"msn" = nc + g*p,"mst" = nc + g*p + g)

ICL = loglik - nc/2*log(nn)

list(ICL=ICL)
}




emmixMOD<- function(p,g,distr,mu,sigma,dof,delta,nrand=10000){
distr<-tolower(distr)
mvrand<-function(n,p,distr,mean,cov,nu,del)
{
switch(distr,
'msn' = rdmsn(n,p,mean=mean,cov=cov,      del=del),
'mst' = rdmst(n,p,mean=mean,cov=cov,nu=nu,del=del),NULL)
}

if(!(distr %in% c("mvn","mvt","msn","mst")))
stop("model specified not available yet")

mu    <- array(mu,c(p,g))
sigma <- array(sigma,c(p,p,g))
dof   <- c(dof)
delta <- array(delta,c(p,g))

modpts <- mu

if(distr %in% c("mvn","mvt"))
return(modpts)


for(h in 1:g) {
dat <-  cbind(mvrand(nrand,p,distr,mu[,h],sigma[,,h],dof[h],delta[,h]))
den <-  ddmix(dat,nrand,p,1, distr,mu[,h],sigma[,,h],dof[h],delta[,h])
id  <-  which.max(den)
modpts[,h] <- dat[id,]
}
return(modpts)
}

emmixfit1<-function(dat,g,clust,distr,ncov,itmax,epsilon,initloop=20)
{
ndist<-switch(tolower(distr),"mvn"=1,"mvt"=2,"msn"=3,"mst"=4,5)
if(ndist>4|ncov<1|ncov>5) 
stop("the model specified is not available yet")

dat<-as.matrix(dat)
n <- nrow(dat)
p <- ncol(dat)

if(n <= 20*g)
stop("sample size is too small")

if(missing(clust) | is.null(clust))
stop("initial clust must be given")

clust <- unclass(as.ordered(clust))

if(max(clust)!=g) stop(paste("The levels of cluster should be g=",g))


obj<-.C('emmixfit1',PACKAGE="JCM",
as.double(dat),as.integer(n),as.integer(p),
as.integer(g),as.integer(ncov),as.integer(ndist),#6
pro   = double(g),mu  = double(p*g),sigma  = double(p*p*g),
dof   = double(g),delta=double(p*g), #11
tau   = double(n*g),double(n*g),double(n*g),double(n*g),double(n*g),
sumtau=double(g), sumvt=double(g),
sumzt=double(g), sumlnv=double(g), 
ewy=double(p*g),ewz=double(p*g),ewyy = double(p*p*g),#23
loglik= double(1),lk= double(itmax),aic= double(1),bic= double(1),
clust = as.integer(clust),#28
error = integer(1),as.integer(itmax),
as.double(epsilon),as.integer(initloop))[c(7:12,24:29)]

lk<-obj$lk;lk<-lk[lk!=0]

list(distr=distr,error=obj$error,loglik=obj$loglik,bic=obj$bic,aic=obj$aic,
pro=obj$pro,mu= array(obj$mu,c(p,g)),
sigma=array(obj$sigma,c(p,p,g)),dof=obj$dof,delta=array(obj$delta,c(p,g)),
clust=obj$clust,tau=array(obj$tau,c(n,g)),lk=lk)
}


emmixfit2<-function(dat,g, init, distr,ncov,itmax,epsilon)
{
ndist<-switch(tolower(distr),"mvn"=1,"mvt"=2,"msn"=3,"mst"=4,5)

if(ndist>4|ncov<1|ncov>5) 
stop("the model specified is not available yet")


dat<-as.matrix(dat)
#
n <- nrow(dat)
p <- ncol(dat)

if(n <= 20*g)
stop("sample size is too small")



if(is.null(init)|missing(init))
stop("init should be provided")

pro   <- init$pro
mu    <- init$mu
sigma <- init$sigma
dof   <- init$dof
delta <- init$delta

obj<-.C('emmixfit2',PACKAGE="JCM",
as.double(dat),as.integer(n),as.integer(p),
as.integer(g),as.integer(ncov),as.integer(ndist),#6
pro   = as.double(pro),mu  = as.double(mu),sigma  = as.double(sigma),
dof   = as.double(dof),delta= as.double(delta), #11
tau   = double(n*g),double(n*g),double(n*g),double(n*g),double(n*g),
sumtau=double(g), sumvt=double(g),
sumzt=double(g), sumlnv=double(g), #20
loglik= double(1),lk= double(itmax),aic= double(1),bic= double(1),
clust = integer(n),#25
error = integer(1),as.integer(itmax),as.double(epsilon))[c(7:12,21:26)]

lk<-obj$lk;lk<-lk[lk!=0]

list(distr=distr,error=obj$error,loglik=obj$loglik,bic=obj$bic,aic=obj$aic,
pro=obj$pro,mu= array(obj$mu,c(p,g)),
sigma=array(obj$sigma,c(p,p,g)),dof=obj$dof,delta=array(obj$delta,c(p,g)),
clust=obj$clust,tau=array(obj$tau,c(n,g)),lk=lk)
}





#*****************************************
#initial values
#*****************************************


initEmmix<-function(dat,g,clust,distr,ncov,maxloop=20)
{
ndist<-switch(tolower(distr),"mvn"=1,"mvt"=2,"msn"=3,"mst"=4,5)

if(ndist>4|ncov<1|ncov>5) 
stop("the model specified is not available yet")

dat<-as.matrix(dat)
n <- nrow(dat)
p <- ncol(dat)

clust<- unclass(as.ordered(clust))


obj<-.C('initfit_',PACKAGE="JCM",
as.double(dat),as.integer(n),as.integer(p),
as.integer(g),as.integer(ncov),as.integer(ndist),
pro   = double(g),mu  = double(p*g),sigma  = double(p*p*g),
dof   = double(g),delta=double(p*g), #11
tau   = double(n*g),double(n*g),double(n*g),double(n*g),double(n*g),
sumtau=double(g), sumvt=double(g),
sumzt=double(g), sumlnv=double(g), 
ewy=double(p*g),ewz=double(p*g),ewyy = double(p*p*g),#23
loglik= double(1),as.integer(clust),
error = integer(1),as.integer(maxloop)) [c(7:11,24,26)]


error <- obj$error
ret   <-NULL

if(error==0) {
ret<-list(distr=distr,error=error,loglik=obj$loglik,
pro=obj$pro,mu= array(obj$mu,c(p,g)),
sigma=array(obj$sigma,c(p,p,g)),
dof=obj$dof,delta=array(obj$delta,c(p,g)))
} else warning("error:",error)

ret

}

init.mix<-function(dat,g,distr,ncov,nkmeans,nrandom,nhclust,maxloop=20)
{
found<-list()
found$loglik<--Inf

n<-nrow(dat)
clust<-rep(1,n)
mclust<-NULL


if(g>1){
if(nkmeans>0) {


for(i in 1:nkmeans)
{
clust<-kmeans(dat,g,nstart=5)$cluster

if(min(table(clust))<10) next

obj<-try(initEmmix(dat,g,clust,distr,ncov,maxloop))

if(length(obj)!=8 | obj$error) next

if(obj$loglik>found$loglik)
{
found<-obj
mclust<-clust
}

}

if(is.null(mclust))
nrandom <- 10
}



if(nrandom>0)
for(i in 1:nrandom) {
clust<-sample(1:g,n,replace=TRUE)

obj<-try(initEmmix(dat,g,clust,distr,ncov,maxloop))
if(length(obj)!=8 | obj$error!=0) next
if(obj$loglik>found$loglik)
{
found<-obj 
mclust<-clust
}

} 

#methods<-c( "ward", "single", "complete", "average", "mcquitty", "median","cen")
methods<-c("complete")
#
if(nhclust)
{
dd <- as.dist((1 - cor(t(dat)))/2)  
#Use correlations between variables ``as distance''

for(j in methods){	
clust<- cutree(hclust(dd, j), k = g)

obj<-try(initEmmix(dat,g,clust,distr,ncov,maxloop))
if(length(obj)!=8 | obj$error!=0) next

if(obj$loglik>found$loglik)
{
found<-obj;
mclust<-clust
}

} #end of for j
} #end if
#------------------------------------------------
}
else
{
obj<-try(initEmmix(dat,g,clust,distr,ncov,maxloop))

if(length(obj)==8)  
{found<-obj;mclust<-clust}

}

if(is.null(mclust)) {
found <-NULL
warning("failed to find a initial values!")
}

found
}





#end


