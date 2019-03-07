

mypanel2<- function(x,y,...) {
par(new=TRUE);
smoothScatter(x,y,..., nrpoints=0)
}

mypanel3<- function(x,y,...) {
par(new=TRUE);
Lab.palette <- colorRampPalette(c("blue", "orange", "red"), space = "Lab")
smoothScatter(x,y, colramp = Lab.palette)

}

mypanel4<- function(x,y,...) {
xy <- cbind(x,y)
par(new=TRUE);
plot(xy, col = densCols(xy), pch=20)
}

panel.density <- function(x, col=1,...)
{
usr <- par("usr"); on.exit(par(usr))
par(usr = c(usr[1:2], 0, 1.5) )

oo = density(x)
y = oo$y
lines(oo$x,y/max(y))

}

conplot <-function(x,y, pro, mu, sigma, dof, delta,distr,
grid=300, nrand=6000,levels=seq(5,95,by=20),col ='white')
{

ndist<-switch(tolower(distr),"mvn"=1,"mvt"=2,"msn"=3,"mst"=4,5)
if(ndist>4) 
stop("the model specified is not available yet")  

ddemmix <- function(dat, n, p, g, distr, pro, mu, sigma, dof, delta)
{
ret<-ddmix(dat,n,p,g, distr, mu,sigma,dof,delta)
c(exp(ret)%*%pro )
} #joint density

    xlim = range(x)+c(-1,0)
    ylim = range(y)+c(-1,0)

g<-length(pro);p<-2

x1 <- seq(xlim[1], xlim[2], length=grid) 
y1 <- seq(ylim[1], ylim[2], length=grid) 

nx <- length(x1)
ny <- length(y1)
xoy <- cbind(rep(x1,ny), as.vector(matrix(y1,nx,ny,byrow=TRUE)))
X <- matrix(xoy, nx*ny, 2, byrow=FALSE)

dens     <- ddemmix(X, nx*ny,2,g, distr, pro, mu, sigma, dof, delta)
dens.mat <- matrix(dens,nx,ny)

n <- table(sample(1:g, nrand, replace = TRUE, prob = pro))
nn <- n
if(length(n) <g) {
nn <- rep(0,g)
for(i in as.numeric(names(n)))
nn[i] <- n[paste(i)]
}
rand     <-  rdemmix(nn,p,g,distr,mu,sigma,dof,delta)
rand.den <-  ddemmix(rand, nrand,2,g, distr, pro, mu, sigma, dof, delta)
cont     <-  quantile(rand.den, prob=levels/100)
contour(x1, y1, dens.mat, levels=cont, add=TRUE,drawlabels=FALSE,lty=1,col =col)
}

conplot2 <- function (x, y, pro, mu, sigma, dof, delta, distr, grid = 300, 
    nrand = 6000, levels = seq(5, 95, by = 20)) 
{
    ndist <- switch(tolower(distr), mvn = 1, mvt = 2, msn = 3, mst = 4, 
        5)
    if (ndist > 4) 
        stop("the model specified is not available yet")
    ddemmix <- function(dat, n, p, g, distr, pro, mu, sigma, 
        dof, delta) {
        ret <- ddmix(dat, n, p, g, distr, mu, sigma, dof, delta)
        c(exp(ret) %*% pro)
    }
    g <- length(pro)
    p <- 2

#make mesh

    xlim = range(x)+c(-1,0)
    ylim = range(y)+c(-1,0)

    x1 <- seq(xlim[1], xlim[2], length= grid)
    y1 <- seq(ylim[1], ylim[2], length= grid)

    nx <- length(x1)
    ny <- length(y1)

    xoy <- cbind(rep(x1, ny), as.vector(matrix(y1, nx, ny, byrow = TRUE)))
    X <- matrix(xoy, nx * ny, 2, byrow = FALSE)
    dens <- ddemmix(X, nx * ny, p, g, distr, pro, mu, sigma, 
        dof, delta)

    dens.mat <- matrix(dens, nx, ny)

#
    n <- table(sample(1:g, nrand, replace = TRUE, prob = pro))
    nn <- n
    if (length(n) < g) {
        nn <- rep(0, g)
        for (i in as.numeric(names(n))) nn[i] <- n[paste(i)]
    }

    rand <- rdemmix(nn, p, g, distr, mu, sigma, dof, delta)
    rand.den <- ddemmix(rand, nrand, 2, g, distr, pro, mu, sigma,dof, delta)
    cont <- quantile(rand.den, prob = 1-levels/100)

    samp <- cbind(x,y)
    samp.den <-ddemmix(samp, length(x), 2, g, distr, pro, mu, sigma,dof, delta)

select <-  which(samp.den>cont)

clust  <-  ifelse(samp.den>cont,2,1)

list(select,clust,x1=x1, y1=y1, dens.mat=dens.mat, cont=cont)
}


conplot3 <- function (x, y, pro, mu, sigma, dof, delta, modpts,distr, grid =300, 
    nrand = 10000, levels = seq(5, 95, by = 20)) 
{
    ndist <- switch(tolower(distr), mvn = 1, mvt = 2, msn = 3, mst = 4, 
        5)
    if (ndist > 4) 
        stop("the model specified is not available yet")
    ddemmix <- function(dat, n, p, g, distr, pro, mu, sigma, 
        dof, delta) {
        ret <- ddmix(dat, n, p, g, distr, mu, sigma, dof, delta)
        c(exp(ret) %*% pro)
    }

    g <- length(pro)
    p <- 2

#----------------------------------------------------
#mesh
#----------------------------------------------------

#make mesh

    xlim = range(x)+c(-1,0)
    ylim = range(y)+c(-1,0)

    x1 <- seq(xlim[1], xlim[2], length=grid)
    y1 <- seq(ylim[1], ylim[2], length=grid)

    nx <- length(x1)
    ny <- length(y1)

    xoy <- cbind(rep(x1, ny), as.vector(matrix(y1, nx, ny, byrow = TRUE)))

    X <- matrix(xoy, nx * ny, 2, byrow = FALSE)

#--------------------------------------------------



for(h in 1:g) { #do each component one by one

    dens <- ddemmix(X, nx * ny, 2, 1, distr, c(1), mu[,h], sigma[,,h],dof[h], delta[,h])

    dens.mat <- matrix(dens, nx, ny)

#  randon sample

    rand <- rdemmix(c(nrand), 2, 1, distr, mu[,h], sigma[,,h], dof[h], delta[,h])

    rand.den <- ddemmix(rand, nrand, 2, 1, distr, c(1), mu[,h], sigma[,,h], dof[h], delta[,h])

#-----------------------------------------

    cont <- quantile(rand.den, prob = 1-levels/100)

    contour(x1, y1, dens.mat, levels = cont, add = TRUE, drawlabels = FALSE,lty = 1, col = h)
    
        if (!is.null(modpts)) 
        points(t(modpts[, h]), col = h,pch=3)
} #end of h loop


}

emmix.filter <- function (S, g=1,distr="mst", diag.panel = TRUE, upper.panel = "type2", 
    lower.panel = "type3", levels = 90, attop = FALSE,title="",path="",plot=TRUE) 
{
    
S <- as.matrix(S)

dat <- S[,1:2]

obj <- emmix(dat,g,distr,ncov=3,itmax=200,debug=0)

ppp <- conplot2(c(dat[,1]), c(dat[,2]), obj$pro, obj$mu, obj$sigma, 
obj$dof, obj$delta, obj$distr, nrand=10000,levels = levels)

clust <- ppp[[2]]
select<- ppp[[1]]

#---------------------------------

mypanel <- function(x, y, ...) {

        par(new = TRUE)
        points(x, y, ..., col = clust)

        st <- pmatch(c(x[1], y[1]), S[1, ])

if(st[1]==1&st[2]==2) {

a=contourLines(ppp$x1, ppp$y1, ppp$dens.mat, levels=ppp$cont)[[1]]

ax <- a$x
#ax[ax<ppp$x1[1]]=ppp$x1[1]+1

ay <- a$y

lines(ax,ay,lty = 2, col = 'blue')

}

}

    if (diag.panel) {
        diag.panel <- panel.density
    }
    else diag.panel <- NULL


    upper.panel <- switch(upper.panel, type2 = mypanel2, type3 = mypanel3, 
        type4 = mypanel4, NULL)

    lower.panel <- switch(lower.panel, type2 = mypanel2, type3 = mypanel3, 
        type4 = mypanel4, NULL)

#-------------------------------

if(plot) {

x11()
pairs(S, pch = ".", panel=mypanel,row1attop = attop, 
        lower.panel = lower.panel, diag.panel = diag.panel,main=paste("Before Filtering (EMMIX):",toupper(distr)))


x11()
pairs(S[select,], upper.panel=upper.panel,row1attop = attop, 
        lower.panel = lower.panel, diag.panel = diag.panel,main=paste("After Filtering (EMMIX):",toupper(distr)))
}

#-------------------------------

if(path!='') {

png(paste(path,'/',title,"-before.png",sep=''),width=512,height=512)
 
pairs(S, pch = ".", panel=mypanel,row1attop = attop, 
        lower.panel = lower.panel, diag.panel = diag.panel,main=paste("Before Filtering (EMMIX):",toupper(distr)))

dev.off()

png(paste(path,'/',title,"-after.png",sep=''),width=512,height=512)
 
pairs(S[select,], upper.panel=upper.panel,row1attop = attop, 
        lower.panel = lower.panel, diag.panel = diag.panel,main=paste("After Filtering (EMMIX):",toupper(distr)))
dev.off()

}

list(subset=select,clust=clust,filter=obj)
}





emmix.contours <- function (S, obj = NULL, clust = NULL,distr="",diag.panel = TRUE, upper.panel = "type2", 
    lower.panel = "type3", levels = seq(5, 95, by = 20), plot=TRUE, title="",path='',attop = FALSE) 
{
    mypanel <- function(x, y, ...) {
        par(new = TRUE)
        smoothScatter(x, y, ..., nrpoints = 0)
        g <- length(obj$pro)
        st <- pmatch(c(x[1], y[1]), S[1, ])
        conplot3(x, y, obj$pro, obj$mu[st, ], obj$sigma[st, st, 
            ], obj$dof, obj$delta[st, ],obj$modpts[st,], obj$distr, levels = levels)
    }
    if (diag.panel) {
        diag.panel <- panel.density
    }
    else diag.panel <- NULL
    upper.panel <- switch(upper.panel, type2 = mypanel2, type3 = mypanel3, 
        type4 = mypanel4, NULL)
    lower.panel <- switch(lower.panel, type2 = mypanel2, type3 = mypanel3, 
        type4 = mypanel4, NULL)


if(plot) {

    if (is.null(clust)) {
        if (!is.null(obj)) {
            pairs(S, panel = mypanel, lower.panel = lower.panel, 
                row1attop = attop, diag.panel = diag.panel,
		main=paste("Contours of Components using EMMIX:",toupper(obj$distr), "Distribution") )
        }
        else {
            pairs(S, upper.panel = upper.panel, lower.panel = lower.panel, 
                row1attop = attop, diag.panel = diag.panel,
		main=paste("The heatmap pairwise plots of the data"))
        }
    }
    else pairs(S, pch = ".", col = clust, row1attop = attop, 
        lower.panel = lower.panel, diag.panel = diag.panel,
	main=paste("Clustering using EMMIX:",toupper(distr), "Distribution") )

}


if(path!='') {

png(paste(path,'/',title,".png",sep=''),width=512,height=512)

    if (is.null(clust)) {
        if (!is.null(obj)) {
            pairs(S, panel = mypanel, lower.panel = lower.panel, 
                row1attop = attop, diag.panel = diag.panel,
		main=paste("Contours of Components using EMMIX:",toupper(obj$distr), "Distribution") )
        }
        else {
            pairs(S, upper.panel = upper.panel, lower.panel = lower.panel, 
                row1attop = attop, diag.panel = diag.panel,
		main=paste("The heatmap pairwise plots of the data"))
        }
    }
    else pairs(S, pch = ".", col = clust, row1attop = attop, 
        lower.panel = lower.panel, diag.panel = diag.panel,
	main=paste("Clustering using EMMIX:",toupper(distr), "Distribution") )

dev.off()

}

}



emmix.flow <-function(S,obj=NULL,distr="",diag.panel=TRUE,
upper.panel="type2",lower.panel="type3",
levels=seq(5,95,by=20),attop=FALSE,clust=NULL,title="",path="",plot=TRUE) {

mypanel<- function(x,y,...) {
par(new=TRUE);
smoothScatter(x,y,..., nrpoints=0)
g <- length(obj$pro)
st <- pmatch(c(x[1],y[1]),S[1,])
conplot(x,y, obj$pro,obj$mu[st,],obj$sigma[st,st,],
obj$dof,obj$delta[st,],obj$distr,levels=levels,nrand=10000)

if(!is.null(obj$modpts))
points(t(obj$modpts[st,]),col=1:g,pch=3)
}



if(diag.panel) {
diag.panel<- panel.density}
else diag.panel<- NULL

upper.panel<-switch(upper.panel,
                    "type2"=mypanel2,
		    "type3"=mypanel3,
		    "type4"=mypanel4,
		    NULL)


lower.panel<-switch(lower.panel,
                    "type2"=mypanel2,
		    "type3"=mypanel3,
		    "type4"=mypanel4,
		    NULL)

if(plot){

    if (is.null(clust)) {
        if (!is.null(obj)) {
            pairs(S, panel = mypanel, lower.panel = lower.panel, 
                row1attop = attop, diag.panel = diag.panel,
		main=paste("Contours of EMMIX:",toupper(obj$distr), "Distribution") )
        }
        else {
            pairs(S, upper.panel = upper.panel, lower.panel = lower.panel, 
                row1attop = attop, diag.panel = diag.panel,
		main=paste("The heatmap pairwise plots of the data"))
        }
    }
    else pairs(S, pch = ".", col = clust, row1attop = attop, 
        lower.panel = lower.panel, diag.panel = diag.panel,
	main=paste("Clustering using EMMIX:",toupper(distr), "Distribution") )



}


if(path!='') {

png(paste(path,'/',title,".png",sep=''),width=512,height=512)

    if (is.null(clust)) {
        if (!is.null(obj)) {
            pairs(S, panel = mypanel, lower.panel = lower.panel, 
                row1attop = attop, diag.panel = diag.panel,
		main=paste("Contours of EMMIX:",toupper(obj$distr), "Distribution") )
        }
        else {
            pairs(S, upper.panel = upper.panel, lower.panel = lower.panel, 
                row1attop = attop, diag.panel = diag.panel,
		main=paste("The heatmap pairwise plots of the data"))
        }
    }
    else pairs(S, pch = ".", col = clust, row1attop = attop, 
        lower.panel = lower.panel, diag.panel = diag.panel,
	main=paste("Clustering using EMMIX:",toupper(distr), "Distribution") )

dev.off()
}

}
