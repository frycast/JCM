
# 1. syntax of unzip data

do.unzip<-function(data,path='.',exdir='.') 
{
#path: A string of the path to where the zipfiles of data are.
#exdir: A string for the name of the project folder where the data are stored.

data <- paste(path,"/",data,".zip",sep='')

unzip(data,exdir=exdir)

}

# 2. syntax to fit individual samples

#fit emmix
do.mix <- function(data="4min",panel="panel1",samples,p=4,g=2,distr="mvt",
init=NULL,itmax=100,header=TRUE,ftype="CSV")
{

#datalist <- dir(paste(data,panel,sep='/'))

#for(datfile in datalist) {

#fname <- tolower(strsplit(toupper(datfile),
#split=paste('.',toupper(ftype),sep=''))[1][1])

for(fname in samples) {

datfile<- paste(fname,ftype,sep='.')

if(tolower(ftype)=="csv")
S <- read.csv(  paste(data,panel,datfile,sep='/'),header=header)[,1:p]

if(tolower(ftype)=="txt")
S <- read.table(paste(data,panel,datfile,sep='/'),header=header)[,1:p]

if(tolower(distr)=="mst") {

if(is.null(init))
{
init       <- emmix(S,g,distr="mvt",itmax=300)
init$dof   <- rep(4,g)
init$delta <- array(0,c(p,g))
}

ret <- try(emmix(S,g,distr=distr,init=init,itmax=itmax))
} 
else
ret <- try(emmix(S,g,distr=distr,itmax=itmax,init=init))

ret$tau <- NULL

if(length(ret)>5) {
dput(ret,paste(data,"_",panel,"_",fname,"_",distr,"_g=",g,".ret",sep=''))

#----------------------
jpeg(paste(data,"_",panel,"_",fname,"_",distr,"_g=",g,".jpeg",sep=''),width=518,height=518)
emmix.contours(S,ret)
dev.off()
}
} #end for

} # end


#fit jcm using previous emmix results as inital values

do.jcm <- function(data="4min",panel="panel1",samples,p=4,g=2,distr="mvt",
itmax=100,header=TRUE,ftype="CSV")
{

#datalist <- dir(paste(data,panel,sep='/'))

#for(datfile in datalist) {

#fname <- tolower(strsplit(toupper(datfile),
#split=paste('.',toupper(ftype),sep=''))[1][1])

for(fname in samples) {

datfile<- paste(fname,ftype,sep='.')

initobj<- dget(paste(data,"_",panel,"_",fname,"_",distr,"_g=",g,".ret",sep=''))

outfile<- paste(data,"_",panel,"_",fname,"_",tolower(distr),
"_g=",g,"_jcm.ret",sep='')

datfile<- paste(data,panel,datfile,sep='/')

initobj$theta<- array(0.1,c(p,g))
initobj$thetu<- rep(0.1,g)

if(tolower(distr)=="mvt") 
ret <- try(jcamvt(datfile, p, g, initobj,itmax=itmax,header=header,ftype=ftype))

if(tolower(distr)=="mst")
ret <- try(jcamst(datfile, p, g, initobj,itmax=itmax,header=header,ftype=ftype))

if(length(ret)>5)
{

if(tolower(ftype)=="csv")
dat <- read.csv(datfile,header=header)[,1:p]
if(tolower(ftype)=="txt")
dat <- read.table(datfile,header=header)[,1:p]

ret$sigma2 <- doentro(p,g, distr, class="jcm", ret) 
ret$modpts <- emmixMOD(p, g, distr, ret$mu, ret$sigma2, 
ret$dof, ret$delta, nrand = 10000)

eobj <- jcamst.estep1(nrow(dat),p,g,dat,ret$pro,ret$mu, 
ret$sigma2, ret$dof, ret$delta)
ret$clust <- tau2clust(eobj$tau)

dput(ret,outfile)

}


} #end for

} # end



# 3. syntax to get distribution templates

do.template <- function(data="4min",panel="panel1",group="LNP_lo",
distr="mvt",p=4,g=6,info,itmax=100,header=TRUE,ftype="csv") 
{

# datafiles

nam <- tolower(info[which(info[,2]==group),1])

datalist <- dir(paste(data,panel,sep='/'))

ids      <- charmatch(nam,datalist)

datfiles <- paste(data,'/',panel,'/',datalist[ids],sep='')

nm       <- length(ids)

# find initials

mm  <- min(nm,10)

files = datfiles[1:mm]

set.seed(12345678)


if(tolower(ftype)=="csv")
dat <- read.csv(files[1],header=header)[,1:p]

if(tolower(ftype)=="txt")
dat <- read.table(files[1],header=header)[,1:p]

n<- nrow(dat)

nn <- sample(1:n,min(n,3000),replace=FALSE)

y  <- dat[nn,]

if(mm>1) {

for(id in 2:mm) {

if(tolower(ftype)=="csv")
dat <- read.csv(files[id],header=header)[,1:p]

if(tolower(ftype)=="txt")
dat <- read.table(files[id],header=header)[,1:p]

n<- nrow(dat)

nn <- sample(1:n,min(n,3000),replace=FALSE)

y  <- rbind(y,dat[nn,])

} #end for

} #end if

if(tolower(distr)=="mst") {

init       <- emmix(y,g,distr="mvt",itmax=300)
init$dof   <- rep(4,g)
init$delta <- array(0,c(p,g))
ret <- try(emmix(y,g,distr=distr,init=init,itmax=itmax))
} 
else
ret <- try(emmix(y,g,distr=distr,itmax=itmax))

dput(ret[1:10],paste("init_",data,"_",panel,"_",group,"_",distr,"_g=",g,".ret",sep=''))

# 1.2 do pooling approach

rm("y")

init <- dget(paste("init_",data,"_",panel,"_",group,"_",distr,"_g=",g,".ret",sep=''))

if(tolower(distr)=="mvt")
obj <- msmvt(datfiles,p,g,init=init,ncov=3,debug=1,itmax=itmax,header=header,ftype=ftype)

if(tolower(distr)=="mst")
obj <- msmst(datfiles,p,g,init=init,ncov=3,debug=1,itmax=itmax,header=header,ftype=ftype)


dput(obj,paste("template_",data,"_",panel,"_",group,"_",distr,"_g=",g,".ret",sep=''))

}

#

JCM.template <- function(data="4min",panel="panel1",group="LNP_lo",
distr="mvt",p=4,g=6,info,itmax=100,header=TRUE,ftype="csv") 
{

# datafiles
nam <- tolower(info[which(info[,2]==group),1])
datalist <- dir(paste(data,panel,sep='/'))
ids      <- charmatch(nam,datalist)
datfiles <- paste(data,'/',panel,'/',datalist[ids],sep='')
# 
init <- dget(paste("template_",data,"_",panel,"_",group,"_",distr,"_g=",g,".ret",sep=''))

init$theta<- array(0.1,c(p,g))
init$thetu<- rep(0.1,g)

if(tolower(distr)=="mvt")
obj <- jcamvt(datfiles,p,g,init=init,ncov=3,debug=1,itmax=itmax,header=header,ftype=ftype)
if(tolower(distr)=="mst")
obj <- jcamst(datfiles,p,g,init=init,ncov=3,debug=1,itmax=itmax,header=header,ftype=ftype)

obj$sigma2 <- doentro(p,g, distr, class="jcm", obj) 
obj$modpts <- emmixMOD(p, g, distr, obj$mu, obj$sigma2, 
obj$dof, obj$delta, nrand = 10000)

dput(obj,paste("template_",data,"_",panel,"_",group,"_",distr,"_g=",g,"_jcm.ret",sep=''))
}

#-------------------------------

do.mould <- function(data,panel,sample,p=4,g=2,
distr="mst",class="") {

dat     <- read.csv(paste(data,'/',panel,'/',sample,".CSV",sep=''))

if(class=="") {
#1) mst

retfile <- paste(data,'_',panel,'_',sample,'_',distr,'_g=',g,".ret",sep='')
obj <- dget(retfile)

mid <- 2
# mould
if( sum( (obj$modpts[2:3,1])^2 ) > sum( (obj$modpts[2:3,2])^2 ))
mid <- 1

mould <- dat[obj$clust == mid,]
}

if(tolower(class)=="jcm") {
#jcm-mst
retfile <- paste(data,'_',panel,'_',sample,'_',distr,'_g=',g,"_jcm.ret",sep='')
jcm <- dget(retfile)

if(is.null(jcm$sigma2)) {

jcm$sigma2 <- doentro(p,g, distr, class="jcm", jcm) 
jcm$modpts <- emmixMOD(p, g, distr, jcm$mu, jcm$sigma2, 
jcm$dof, jcm$delta, nrand = 10000)

eobj <- jcamst.estep1(nrow(dat),p,g,dat,jcm$pro,jcm$mu, jcm$sigma2, jcm$dof, jcm$delta)
jcm$clust <- tau2clust(eobj$tau)
}

mid <- 2
# mould
if( sum( (jcm$modpts[2:3,1])^2 ) > sum( (jcm$modpts[2:3,2])^2 ))
mid <- 1

mould <- dat[jcm$clust == mid,]

}

mould

}

#-------------------------------------------


doentro<- function(p,g, distr, class, obj) 
{
distr<- tolower(distr)

if(is.na(match(distr,c("mvn","mvt","msn","mst"))))
stop("distribution specified is not available!")

    mu <- array(obj$mu, c(p, g))
    sigma <- array(obj$sigma, c(p, p, g))
    delta <- array(0, c(p, g))
    dof <- rep(0, g)
    if (distr == "mst") {
        delta <- array(obj$delta, c(p, g))
        dof <- obj$dof
    }
    if (distr == "mvt") {
        dof <- obj$dof
    }
    if (distr == "msn") {
        delta <- array(obj$delta, c(p, g))
    }
    if (tolower(class) == "jcm") {
        theta <- array(obj$theta, c(p, g))
        thetu <- obj$thetu
        omiga <- array(0, c(p, p, g))
        onemat <- array(1, c(p, p))
        for (h in 1:g) {
            z1 <- diag(mu[, h], p)
            omiga[, , h] <- (sigma[, , h] + array((z1 * theta[, 
                h]) %*% t(z1), c(p, p)) + thetu[h] * onemat)
        }
        sigma <- omiga
    }
sigma
}

#----------------
