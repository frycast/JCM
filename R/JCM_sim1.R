
# setting 1: changes of number of observations N

# data folder: JCM/jcm-sim/simulation1/n010, n020,n030,...,n100.

# sample names: simdat-01.txt,02,...,15 samples on each setting

#setwd("JCM")  #first

#system("mkdir jcm-sim")

#setwd("jcm-sim") #second

#system("mkdir simulation1")
#system("mkdir simulation2")
#system("mkdir simulation3")
#system("mkdir simulation4")

#simulation 1
#setwd("simulation1")

#system("mkdir n010")
#system("mkdir n020")
#system("mkdir n030")
#system("mkdir n040")
#system("mkdir n050")
#system("mkdir n060")
#system("mkdir n070")
#system("mkdir n080")
#system("mkdir n090")
#system("mkdir n100")

#setwd("../../../")


#setwd("d:/uq/uq2011")
#source("EmSkew.R")

dosim1<-function(nk=10,path="JCM/jcm-sim/simulation1")
{

datpath =  paste(path,"/n0",nk,sep='')

if(nk==100)
datpath =  paste(path,"/n100",sep='')

p=4
g=4

n=nk*1000

#--------------------------------------

# set parameters

mu  <- cbind(c(0,0,5,5),c(5,5,0,0),c(0,0,0,0),c(5,0,5,0))
sigma <- array(0,c(p,p,g))
sigma[,,1] <- diag(4)*0.5
sigma[,,2] <- diag(4)*0.5
sigma[,,3] <- diag(4)*0.5
sigma[,,4] <- diag(4)*0.5

#dof <- c(4,5,6,7)

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

for(j in 1:10) {


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

delta <- c(0.1,0.15,0.2,0.23,0.25)

for(j in 1:5) {

pro <- c(0.25-delta[j],0.25+delta[j],rep(0.25,2))


nn <- table(sample(1:g, n, replace = TRUE, prob = pro))
names(nn) <- NULL

if(length(nn)<g)
nn <- c(0,nn)

y   <- rdemmix(nn,p,g,"mvn",mue,sigma)

ra  <- rdmvn(n,p,mua,theta)

rb  <- rdemmix(nnb,p,g,"mvn",mub,thetu)

count=0
for(h in 1:g) {
if(nn[h]>0) {
ra[count+1:nn[h],]<- t(t(ra[count+1:nn[h],])*mu[,h]+rb[h,])
count=count+nn[h]
}}

x  <-  y+ra

#x11();pairs(x,col=rep(1:g,nn))}


datfile <- paste(datpath,"/simdat-",10+j,".txt",sep='')
write.table(x,datfile,row.names=FALSE,col.names=paste("v",1:p,sep='') )

} # end of loop

#-----------------------------


}



dojcm1 <- function(nk=10,dpath="JCM/jcm-sim/simulation1",
opath="JCM/jcm-sim")
{

datpath =  paste(dpath,"/n0",nk,sep='')

if(nk==100)
datpath =  paste(dpath,"/n100",sep='')

t1 <- proc.time()


p=4
g=4

n=nk*1000

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

ret <- jcamst(datfiles,p,g,init,ncov=3,dofon=0,
debug=1,itmax=100,epsilon=1e-10,ftype='txt' )

t2 <- proc.time()-t1

ret$time <- t2

dput(ret,paste(opath,"/simulation1-",nk,".ret",sep=''))

}





