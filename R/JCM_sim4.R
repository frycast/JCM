
# setting 4: changes of number of cohort samples from 5 to 50 by 5

# data folder: JCM/jcm-sim/simulation4/ks05, ks10,ks15,...,ks50.

# sample names: simdat-01.txt,02,...,15 samples on each setting

#setwd("JCM/jcm-sim")  #first

#system("mkdir simulation4")

#simulation 4
#setwd("simulation4")

#system("mkdir ks05")
#system("mkdir ks10")
#system("mkdir ks15")
#system("mkdir ks20")
#system("mkdir ks25")
#system("mkdir ks30")
#system("mkdir ks35")
#system("mkdir ks40")
#system("mkdir ks45")
#system("mkdir ks50")

#setwd("../../../")

#setwd("d:/uq/uq2011")
#source("EmSkew.R")



dosim4<-function(ks=10,path="JCM/jcm-sim/simulation4")
{

if(ks<10)
datpath =  paste(path,"/ks0",ks,sep='')

if(ks>=10)
datpath =  paste(path,"/ks",ks,sep='')

p=4
g=4
n=20000

#--------------------------------------

# set parameters

mu  <- cbind(c(0,0,5,5),c(5,5,0,0),c(0,0,0,0),c(5,0,5,0))
sigma <- array(0,c(p,p,g))
sigma[,,1] <- diag(4)*0.5
sigma[,,2] <- diag(4)*0.5
sigma[,,3] <- diag(4)*0.5
sigma[,,4] <- diag(4)*0.5

# theta_a : random a at cell-level

mua  = rep(1,p)
theta= diag(p)*0.1

#ra  <- rdmvn(n, p,mua,theta)

# theta_b at cluster level b

thetu <- array(0,c(p,p,g))
mub = array(0,c(p,g))
thetu[,,1] = diag(p)*0.1 
thetu[,,2] = diag(p)*0.15
thetu[,,3] = diag(p)*0.2
thetu[,,4] = diag(p)*0.25

nnb = rep(1,g)

#rb  <- rdemmix(nnb,p,g,"mvn",mub,thetu)

#-------------------
mue <- (array(0,c(p,g)))
pro <- rep(0.25,4)

nn <- table(sample(1:g, n, replace = TRUE, prob = pro))
names(nn) <- NULL

for(j in 1:ks) {


y   <- rdemmix(nn,p,g,"mvn",mue,sigma)

ra  <- rdmvn(n, p,mua,theta)

rb  <- rdemmix(nnb,p,g,"mvn",mub,thetu)

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

}



dojcm4 <- function(ks=10,dpath="JCM/jcm-sim/simulation4",
opath="JCM/jcm-sim")
{

t1 <- proc.time()

if(ks<10)
datpath =  paste(dpath,"/ks0",ks,sep='')

if(ks>=10)
datpath =  paste(dpath,"/ks" ,ks,sep='')

p=4
g=4
n=20000

# set initials

mu  <- cbind(c(0,0,5,5),c(5,5,0,0),c(0,0,0,0),c(5,0,5,0))
sigma <- array(0,c(p,p,g))
sigma[,,1] <- diag(4)*0.5
sigma[,,2] <- diag(4)*0.5
sigma[,,3] <- diag(4)*0.5
sigma[,,4] <- diag(4)*0.5

dof <- c(4,4,4,4)
delta <- array(0,c(p,g))

theta= array(0.1,c(p,g))
thetu <- c(0.1,0.15,0.2,0.25)

init       <- list()

init$mu    <- mu
init$sigma <- sigma
init$dof   <- dof
init$delta <- delta
init$theta <- theta
init$thetu <- thetu
init$pro   <- rep(0.25,g)
#-------------------------

datfiles <- dir(datpath,full.names=TRUE)

ret <- jcamst(datfiles,p,g,init,ncov=3,dofon=0,debug=1,
itmax=100,epsilon=1e-10,ftype='txt' )

t2 <- proc.time()-t1

ret$time <- t2

dput(ret,paste(opath,"/simulation4-",ks,".ret",sep=''))

}





