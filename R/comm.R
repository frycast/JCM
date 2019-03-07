rand.index <- function (x, y) 
{
#This function is a modified copy of mclust::adjustedRandIndex
#pls invoke mclust pacakge to see the copy right and refer to mclust manual.pdf.

x <- unclass(as.ordered(x))
y <- unclass(as.ordered(y))

    
    x <- as.vector(x)
    y <- as.vector(y)
    xx <- outer(x, x, "==")
    yy <- outer(y, y, "==")
    upper <- row(xx) < col(xx)
    xx <- xx[upper]
    yy <- yy[upper]
    a <- sum(as.numeric(xx & yy))
    b <- sum(as.numeric(xx & !yy))
    c <- sum(as.numeric(!xx & yy))
    d <- sum(as.numeric(!xx & !yy))
    ni <- (b + a)
    nj <- (c + a)
    abcd <- a + b + c + d
    q <- (ni * nj)/abcd
#  rand index
    rd1<-(a+d)/abcd
#  adjusted rand index
    rda<-(a - q)/((ni + nj)/2 - q)
ind <- c(rd1,rda)
names(ind) <- c("Rand Index","Adjusted Rand Index")
ind
}

#------------------


inverse <-function(sigma,p)
{
if(length(c(sigma))!=p*p | ncol(sigma)!=p)
stop("sigma should be p by p matrix")

obj <- .Fortran('inverse3',PACKAGE="JCM",
as.double(sigma),inv=double(p*p), 
det=double(1),as.integer(p), error = integer(1),
count = integer(1),index = integer(p))

if(obj$error) stop("")
a <- array(obj$inv,c(p,p))
as.matrix(t(a)%*%a)
}
tau2clust<-function(tao)
{
apply(tao,FUN=which.max,MARGIN=1)
}

getcov <-function(msigma,sumtau,n,p,g,ncov)
{
sigma<-array(0,c(p,p))
if( (ncov==1)|(ncov==2))
{
for(h in 1:g)
sigma<-sigma+sumtau[h]*msigma[,,h]
sigma<-as.matrix(sigma/n)

if(ncov==2)
sigma<-diag(c(diag(sigma)),p)
for(h in 1:g)
msigma[,,h]=sigma
}

if(p>1)
{
if(ncov==4)
for(h in 1:g)
msigma[,,h]<-diag(c(diag(msigma[,,h])),p)

if(ncov==5)
for(h in 1:g)
msigma[,,h]<-diag(sum(diag(msigma[,,h]))/p,p)
}

msigma
}


mvt.dof <-
function(sumtau,sumlnv,lx=2+1e-4,ux=200)
{

if(sumtau <=2)
return(4L)

f<-function(v,sumlnv,sumtau) 
{
sumtau*( log(v/2)-digamma(v/2)+1)+ sumlnv
}

if(f(lx,sumlnv,sumtau)*f(ux,sumlnv,sumtau)>0)
return(ux)
else
(uniroot(f,c(lx,ux),sumlnv=sumlnv,sumtau=sumtau)$root)
}



error.rate<-function(clust1,clust2)
{

clust1 <- unclass(as.ordered(clust1))
clust2 <- unclass(as.ordered(clust2))

if((n=length(clust1))!=length(clust2))
{warning("error: length not equal");return}
if( (g=length(table(clust1)))!=length(table(clust2)))
{warning("the number of clusters are not equal");return}

permute<-function(a)
{
n<-length(a)
if(n==1)
f<-a
else
{
nm<-gamma(n)
f<-array(0,c(n,n*nm))
j<-1

for(i in a)
{
 f[1, (j-1)*nm+1:nm]<-i
 f[-1,(j-1)*nm+1:nm]<-permute(setdiff(a,i))
 j<-j+1
}
}

f
}


#
id<-1:n
cmb<-permute(1:g)
nperm<-ncol(cmb)
rate<-rep(0,nperm)
#
for(i in 1:nperm)
{
tmp<-rep(0,g)
tc<-rep(0,n)
for(j in 1:g)
tc[clust2==j]=cmb[j,i]

for(j in 1:g)
{  
tmp1<-0 
for(k in (1:g)[-j])
        tmp1<-tmp1+length(intersect(id[clust1==j],id[tc==k]))
tmp[j]<-tmp1
}
rate[i]<-sum(tmp)/n
}
min(rate)
}


#--------------------------------------------

# 0's have to be removed in this version

F.Measures<- function(Clabels,Klabels) {

Clab <- unclass(as.ordered(Clabels))
Klab <- unclass(as.ordered(Klabels))

if((n=length(Clab))!=length(Klab))
{warning("does not match");return}

g=max(Clab)  #true number of clusters
k=max(Klab)  #user defined number of clusters

nc <- rep(0,g)
nk <- rep(0,k)

for(i in 1:g) 
nc[i] <- sum(Clab==i)

for(j in 1:k) 
nk[j] <- sum(Klab==j)

#------------------------
fvalue<-0

for(i in 1:g) {
C.index <- which(Clab==i)
a<- rep(0,k)
for(j in 1:k) {
K.index <- which(Klab==j)
a[j] <- 2*sum(C.index %in% K.index)/(nc[i]+nk[j])
}
fvalue<- fvalue+nc[i]*max(a)/n
}
fvalue
}


F.CI <- function(fvalues,B=1e4){
vmean <- rep(0,B)
nk <- length(fvalues)
for(it in 1:B) {
x <- sample(fvalues,nk,replace=TRUE)
vmean[it] <- mean(x)
}
list(mean=mean(vmean), CI95=quantile(vmean,c(0.05,0.95)))
}

#example:
#fvalues <- c(0.88,0.98,0.88,0.56,0.55,0.23,0.44,0.55,0.77,0.77,0.88)
#F.CI(fvalues)

#----------------------------

#--------------------------------------

# the end
