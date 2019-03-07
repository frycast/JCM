
# setting 3: changes of the number of components g from 1 to 10 by 1

# data folder: JCM/jcm-sim/simulation3/g01,g02,...,g10.

# sample names: simdat-01.txt,02,...,15 samples on each setting

#setwd("d:/uq/uq2011/JCM/jcm-sim")

#system("mkdir simulation3")

#simulation 3
#setwd("simulation3")

#system("mkdir g01")
#system("mkdir g02")
#system("mkdir g03")
#system("mkdir g04")
#system("mkdir g05")
#system("mkdir g06")
#system("mkdir g07")
#system("mkdir g08")
#system("mkdir g09")
#system("mkdir g10")

#setwd("d:/uq/uq2011")
#source("EmSkew.R")

dosim3<-function(ng=1,path="JCM/jcm-sim/simulation3")
{

if(ng<10) 
datpath =  paste(path,"/g0",ng,sep='')

if(ng>=10)
datpath =  paste(path,"/g",ng,sep='')

p=4
g=ng
n=20000

#--------------------------------------

# set parameters

if(g==1)
mu = cbind(c(0,0,0,0))


if(g==2)
mu = cbind(c(0,0,0,0),c(6,6,0,0))


if(g==3) 
mu = cbind(c(0,0,0,0),c(0,6,0,0),c(6,0,6,0))


if(g==4)
mu = cbind(c(0,0,0,0),c(0,6,0,0),
c(6,0,0,0),c(6,6,0,0))

if(g==5)
mu = cbind(c(0,0,0,0),
mu2 = c(0,6,0,0),
mu3 = c(6,0,0,0),
mu4 = c(6,6,0,0),
mu5 = c(0,0,5,0))

if(g==6)
mu =cbind( c(0,0,0,0),
mu2 = c(0,6,0,0),
mu3 = c(6,0,0,0),
mu4 = c(6,6,0,0),
mu5 = c(0,0,5,0),
mu6 = c(6,0,5,0))

if(g==7)
mu = cbind(c(0,0,0,0),
mu2 = c(0,6,0,0),
mu3 = c(6,0,0,0),
mu4 = c(6,6,0,0),
mu5 = c(0,0,5,0),
mu6 = c(6,0,5,0),
mu7 = c(6,6,5,0))

if(g==8)
mu = cbind(c(0,0,0,0),
mu2 = c(0,6,0,0),
mu3 = c(6,0,0,0),
mu4 = c(6,6,0,0),
mu5 = c(0,0,5,0),
mu6 = c(6,0,5,0),
mu7 = c(6,6,5,0),
mu8 = c(0,6,5,0))

if(g==9)
mu = cbind(c(0,0,0,0),
mu2 = c(0,6,0,0),
mu3 = c(6,0,0,0),
mu4 = c(6,6,0,0),
mu5 = c(0,0,5,0),
mu6 = c(6,0,5,0),
mu7 = c(6,6,5,0),
mu8 = c(0,6,5,0),
mu9 = c(0,0,0,5))

if(g==10)
mu = cbind(c(0,0,0,0),
mu2 = c(0,6,0,0),
mu3 = c(6,0,0,0),
mu4 = c(6,6,0,0),
mu5 = c(0,0,5,0),
mu6 = c(6,0,5,0),
mu7 = c(6,6,5,0),
mu8 = c(0,6,5,0),
mu9 = c(0,0,0,5),
mu10 = c(0,0,5,5))

#--------------------------------------

sigma <- array(0,c(p,p,g))

for(h in 1:g) 
sigma[,,h] <- diag(p)*0.5

# theta_a : random a at cell-level

mua  = rep(1,p)
theta= diag(p)*0.1

#ra  <- rdmvn(n, p,mua,theta)

# theta_b at cluster level b

thetu <- array(0,c(p,p,g))
mub = array(0,c(p,g))

for(h in 1:g) 
thetu[,,h] = diag(p)*0.2 

nnb = rep(1,g)

#rb  <- rdemmix(nnb,p,g,"mvn",mub,thetu)

#-------------------
mue <- (array(0,c(p,g)))
pro <- rep(1/g,g)

nn <- table(sample(1:g, n, replace = TRUE, prob = pro))
names(nn) <- NULL

if((dd=sum(nn)-n)!=0)
nn[1]=nn[1]-dd

for(j in 1:15) {

y   <- cbind(rdemmix(nn,p,g,"mvn",mue,sigma))

ra  <- cbind(rdmvn(n, p,mua,theta))

rb  <- array(rdemmix(nnb,p,g,"mvn",mub,thetu),c(g,p))

count=0
for(h in 1:g) {
ra[count+1:nn[h],]<- t(t(ra[count+1:nn[h],])*mu[,h]+rb[h,])
count=count+nn[h]
}

x  <-  y+ra

#x11();pairs(x,col=rep(1:g,nn))}

datfile <- paste(datpath,"/simdat-",j,".txt",sep='')
write.table(x,datfile,row.names=FALSE,col.names=paste("v",1:p,sep='') )

} # end of loop


#-----------------------------


}



dojcm3 <- function(ng=2,dpath="JCM/jcm-sim/simulation3",
opath="JCM/jcm-sim") {

t1 <- proc.time()

# set initials

#-------------------------------

if(ng<10) 
datpath =  paste(dpath,"/g0",ng,sep='')

if(ng>=10)
datpath =  paste(dpath,"/g",ng,sep='')

p=4
g=ng
n=20000

#--------------------------------------

# set parameters

if(g==1)
mu = cbind(c(0,0,0,0))


if(g==2)
mu = cbind(c(0,0,0,0),c(6,6,0,0))


if(g==3) 
mu = cbind(c(0,0,0,0),c(0,6,0,0),c(6,0,6,0))


if(g==4)
mu = cbind(c(0,0,0,0),c(0,6,0,0),
c(6,0,0,0),c(6,6,0,0))

if(g==5)
mu = cbind(c(0,0,0,0),
mu2 = c(0,6,0,0),
mu3 = c(6,0,0,0),
mu4 = c(6,6,0,0),
mu5 = c(0,0,5,0))

if(g==6)
mu =cbind( c(0,0,0,0),
mu2 = c(0,6,0,0),
mu3 = c(6,0,0,0),
mu4 = c(6,6,0,0),
mu5 = c(0,0,5,0),
mu6 = c(6,0,5,0))

if(g==7)
mu = cbind(c(0,0,0,0),
mu2 = c(0,6,0,0),
mu3 = c(6,0,0,0),
mu4 = c(6,6,0,0),
mu5 = c(0,0,5,0),
mu6 = c(6,0,5,0),
mu7 = c(6,6,5,0))

if(g==8)
mu = cbind(c(0,0,0,0),
mu2 = c(0,6,0,0),
mu3 = c(6,0,0,0),
mu4 = c(6,6,0,0),
mu5 = c(0,0,5,0),
mu6 = c(6,0,5,0),
mu7 = c(6,6,5,0),
mu8 = c(0,6,5,0))

if(g==9)
mu = cbind(c(0,0,0,0),
mu2 = c(0,6,0,0),
mu3 = c(6,0,0,0),
mu4 = c(6,6,0,0),
mu5 = c(0,0,5,0),
mu6 = c(6,0,5,0),
mu7 = c(6,6,5,0),
mu8 = c(0,6,5,0),
mu9 = c(0,0,0,5))

if(g==10)
mu = cbind(c(0,0,0,0),
mu2 = c(0,6,0,0),
mu3 = c(6,0,0,0),
mu4 = c(6,6,0,0),
mu5 = c(0,0,5,0),
mu6 = c(6,0,5,0),
mu7 = c(6,6,5,0),
mu8 = c(0,6,5,0),
mu9 = c(0,0,0,5),
mu10 = c(0,0,5,5))

#--------------------------------------

sigma <- array(0,c(p,p,g))

for(h in 1:g) 
sigma[,,h] <- diag(p)*0.5

#-------------------------------

dof <- rep(4,g)
delta <- array(0,c(p,g))

theta= array(0.1,c(p,g))
thetu <- rep(0.2,g)

init       <- list()
init$mu  <- mu
init$sigma <- sigma
init$dof   <- dof
init$delta <- delta
init$theta <- theta
init$thetu <- thetu
init$pro <- rep(0.25,g)
#-------------------------

datfiles <- dir(datpath,full.names=TRUE)

ret <- jcamst(datfiles,p,g,init,ncov=3,dofon=0,
debug=1,itmax=100,epsilon=1e-10,ftype='txt')

t2 <- proc.time()-t1

ret$time <- t2

dput(ret,paste(opath,"/simulation3-",ng,".ret",sep=''))

}





