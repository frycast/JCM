do.kld <-function (datfile, SelfTemplate, ClassTemplates, st, distr = "mvt", 
    class = "jcm", header = TRUE,ftype="csv") 
{
    s2cReKLD <- function(datfile, distr, class, retfile, st) {

if(ftype=="txt")
        dat <- read.table(datfile, header = header)[, st]
if(ftype=="csv")
        dat <- read.csv(datfile, header = header)[, st]

obj <- dget(retfile)
entropy(dat,distr,class,obj) 
}
d=length(ClassTemplates)
kld <- rep(0,d)
kld0   <- s2cReKLD(datfile,distr,class,SelfTemplate,     st)
for(j in 1:d) {
kld[j] <- s2cReKLD(datfile,distr,class,ClassTemplates[j],st)
}
kld-kld0
}



entropy <-function(dat,distr,class,obj) {

S <- as.matrix(dat)
p <- ncol(S)
n <- nrow(S)

# format parameters
pro <- obj$pro
g   <- length(pro)
mu    <- array(obj$mu   ,c(p,g)  )
sigma <- array(obj$sigma,c(p,p,g))
delta <- array(0,c(p,g))
dof   <- rep(0,g)

if(tolower(distr)=="mst") {
delta <- array(obj$delta,c(p,g))
dof   <- obj$dof
}
if(tolower(distr)=="mvt") {
dof<- obj$dof}
if(tolower(distr)=="msn") {
delta <- array(obj$delta,c(p,g))}

# adjust JCM sigma
if(class=="jcm") {
theta <- array(obj$theta,c(p,g)  )
thetu <-       obj$thetu
omiga <- array(0,c(p,p,g))
onemat <- array(1,c(p,p))  #z2 %*% t(z2)
for(h in 1:g) {
z1 <- diag(mu[,h],p)
omiga[,,h]<-(sigma[,,h]+array((z1*theta[,h])%*%t(z1),c(p,p))+thetu[h]*onemat)
}
sigma <- omiga
}

# likelihood
logden  <- ddmix(S,n,p,g,distr, mu,sigma,dof,delta)
(-mean(log(c(exp(logden)%*%pro ))))
}
