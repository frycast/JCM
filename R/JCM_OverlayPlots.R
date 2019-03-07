

# NOTICE

# The following program is a modified verion
# of program "plotmixt.skew.6.R" from Pyne 's team at MIT.
# to fit into emmix.
#

###############################################################################
## Creates plots of mixture density functions
#
## Parameters
## mu - means
## sigma - variances
## pro - vector of proportions of each mixture component 
## dof - degrees of freedom
## distr - "mvt" - normal mixture
##      - "mvt" - t mixture; "mst" - skew t-distribution;
##      "msn" - skew normal distribution
## ...
###############################################################################


plotmixt.skew = function(retfiles,which.dim,distr="mvt", 
cont=c(25,50,75),  xlim, ylim, zlim, gridsize, 
nrand=1e5, var.label, line.col, marginal.col,  
feature.col, lwd=1,...)
{
add.feature=TRUE

marginal.scale=1
marginal.ann.col="blue"
supp=3.7


###############################################################################
# Permute a list of values
#
# Same function as EXPAND.GRID (base package), modified to take 
# list as an argument and returns a matrix 
###############################################################################



permute <- function(args) 
{

nargs <- length(args)
if (!nargs) return(as.data.frame(list()))

a1 <- args[[1]]

if (nargs == 1 && is.list(a1)) {
	nargs <- length(a1)
	args <- a1
}

  if (nargs <= 1) 
    return(as.data.frame(if (nargs == 0 || is.null(args[[1]])) list() else args, 
                         optional = TRUE))
  cargs <- args
  rep.fac <- 1
  ##orep = final.len = prod(sapply(args, length))
  orep <- prod(sapply(args, length))
  
  for (i in 1:nargs) {
    x = args[[i]]
    nx = length(x)
    orep = orep/nx
    cargs[[i]] = rep(rep(x, rep(rep.fac, nx)), orep)
    rep.fac = rep.fac * nx
  }
  do.call("cbind", cargs)
} 

convert.dget <- function(file, which.dim)
{
  pm.list = list() 
  j = 1
  for (f in file)
  {  
    params = dget(f)
    
    mu = (params$mu[which.dim,])

    sigma = params$sigma[which.dim,which.dim,]

    dof = params$dof

    pro = round(params$pro, 3)

    delta = params$delta[which.dim,]

    pm.list[[j]] = list(mu=mu, sigma=sigma, dof=dof, pro=pro, delta=delta)
    j = j+1
  }

  return(pm.list)
}

ddemmix <- function(dat, n, p, g, distr, pro, mu, sigma, dof, delta)
{
ret<-ddmix(dat,n,p,g, distr, mu,sigma,dof,delta)
c(exp(ret)%*%pro )
} #joint density

params<- convert.dget(retfiles, which.dim)


  if (missing(line.col)) line.col = 1:length(params)
  if (missing(marginal.col)) marginal.col = line.col
  if (missing(feature.col)) feature.col = rep(1, length(line.col))
 
  j = 0

  dens.list = numeric()
  dens.cont.list = numeric()
  fs.list = numeric()
  
  for (pm in params)
  {
    mu    = pm$mu
    pro   = pm$pro
    sigma = pm$sigma
    delta = pm$delta
    dof   = pm$dof
    j = j+1

g <- length(pro)
d = 3


aaa<-c()
for(h in 1:g) 
aaa <- c(aaa,sqrt(diag(sigma[,,h])))




    maxSigmas = supp*max(aaa)
    
    w = rep(1, ncol(mu))

    if (missing(var.label)) var.label = c("x", "y", "z")

    if (missing(xlim)) xlim = c(min(mu[1,]) - maxSigmas, max(mu[1,]) + maxSigmas)
    if (missing(ylim)) ylim = c(min(mu[2,]) - maxSigmas, max(mu[2,]) + maxSigmas)
    if (missing(zlim)) zlim = c(min(mu[3,]) - maxSigmas, max(mu[3,]) + maxSigmas)

    if (missing(gridsize)) gridsize = rep(51,d)
    if (missing(line.col)) line.col = 1 

    x1 = seq(xlim[1], xlim[2], length=401)
    y1 = seq(ylim[1], ylim[2], length=401)
    z1 = seq(zlim[1], zlim[2], length=401)

    x = seq(xlim[1], xlim[2], length=gridsize[1])
    y = seq(ylim[1], ylim[2], length=gridsize[2])
    z = seq(zlim[1], zlim[2], length=gridsize[3])

    xy = permute(list(x,y))
    xz = permute(list(x,z))
    yz = permute(list(y,z))

	nx <- length(x1)
	ny <- length(y1)
	nz <- length(z1)



	xdens = ddemmix(x1,nx,1,g,distr,pro,mu[1,],sigma[1,1,],dof,delta[1,])
	x.rand = rdemmix2(10*nrand,1,g,distr,pro,mu[1,],sigma[1,1,],dof,delta[1,])
	if (add.feature)
	{
	xh = hmise.mixt(samp=10*nrand, mus=0, sigma=sd(x.rand), props=1, deriv.order=2)
	xfs = featureSignif(x=x.rand, bw=xh)
	x=xfs$fhat$x.grid[[1]]
	xdens2 =ddemmix(x,length(x),1,g,distr,pro,mu[1,],sigma[1,1,],dof,delta[1,])
	xfhat = list(x.grid=xfs$fhat$x.grid, est=xdens2)
	}


	ydens = ddemmix(y1,ny,1,g,distr,pro,mu[2,],sigma[2,2,],dof,delta[2,])
	y.rand = rdemmix2(10*nrand,1,g,distr,pro,mu[2,],sigma[2,2,],dof,delta[2,])
	if (add.feature)
	{
	yh = hmise.mixt(samp=10*nrand, mus=0, sigma=sd(y.rand), props=1, deriv.order=2)
	yfs = featureSignif(x=y.rand, bw=yh)
	x=yfs$fhat$x.grid[[1]]
	ydens2 =ddemmix(x,length(x),1,g,distr,pro,mu[2,],sigma[2,2,],dof,delta[2,])
	yfhat = list(x.grid=yfs$fhat$x.grid, est=ydens2) 
	}

	zdens = ddemmix(z1,nz,1,g,distr,pro,mu[3,],sigma[3,3,],dof,delta[3,])
	z.rand = rdemmix2(10*nrand,1,g,distr,pro,mu[3,],sigma[3,3,],dof,delta[3,])
	if (add.feature)
	{
	zh = hmise.mixt(samp=10*nrand, mus=0, sigma=sd(z.rand), props=1, deriv.order=2)
	zfs = featureSignif(x=z.rand, bw=zh)
	x=zfs$fhat$x.grid[[1]]
	zdens2 =ddemmix(x,length(x),1,g,distr,pro,mu[3,],sigma[3,3,],dof,delta[3,])
	zfhat = list(x.grid=zfs$fhat$x.grid, est=zdens2) 
	}


#--------------------------

st=c(1,2)

	xydens = ddemmix(xy,nrow(xy),2,g,distr,pro,mu[st,],sigma[st,st,],dof,delta[st,])
	xy.rand = rdemmix2(nrand,2,g,distr,pro,mu[st,],sigma[st,st,],dof,delta[st,])
	xydens.rand = ddemmix(xy.rand,nrand,2,g,distr,pro,mu[st,],sigma[st,st,],dof,delta[st,])

st=c(1,3)

	xzdens = ddemmix(xz,nrow(xy),2,g,distr,pro,mu[st,],sigma[st,st,],dof,delta[st,])
	xz.rand = rdemmix2(nrand,2,g,distr,pro,mu[st,],sigma[st,st,],dof,delta[st,])
	xzdens.rand = ddemmix(xz.rand,nrand,2,g,distr,pro,mu[st,],sigma[st,st,],dof,delta[st,])
  
st=c(2,3)

	yzdens = ddemmix(yz,nrow(yz),2,g,distr,pro,mu[st,],sigma[st,st,],dof,delta[st,])
	yz.rand = rdemmix2(nrand,2,g,distr,pro,mu[st,],sigma[st,st,],dof,delta[st,])
	yzdens.rand = ddemmix(yz.rand,nrand,2,g,distr,pro,mu[st,],sigma[st,st,],dof,delta[st,])

#--------------------------


    dens = rbind(cbind(xz, xzdens, 1), cbind(yz, yzdens, 2), cbind(z1, z1, zdens, 3), cbind(xy, xydens, 4), cbind(y1, y1, ydens, 5),  cbind(yz[,2:1], yzdens, 6), cbind(x1, x1, xdens,7), cbind(xy[,2:1], xydens, 8), cbind(xz[,2:1], xzdens, 9))
    dens = data.frame(x=dens[,1], y=dens[,2], f=dens[,3], lab=dens[,4])
    dens$lab = as.factor(dens$lab)
    dens.list = rbind(dens.list, cbind(dens, which.mixt=j))

    xyhts = quantile(xydens.rand, prob=(100-cont)/100)
    xzhts = quantile(xzdens.rand, prob=(100-cont)/100)
    yzhts = quantile(yzdens.rand, prob=(100-cont)/100)

    dens.cont = rbind(cbind(xzhts, 1), cbind(yzhts, 2), cbind(xyhts, 4), cbind(yzhts, 6), cbind(xyhts, 8), cbind(xzhts, 9))
    dens.cont = data.frame(cont=dens.cont[,1], lab=dens.cont[,2])
    dens.cont$lab = as.factor(dens.cont$lab)
    dens.cont.list = rbind(dens.cont.list, cbind(dens.cont, which.mixt=j))

#-------------------------


    if (add.feature)
    {
      fs = rbind(cbind(zfhat$x.grid[[1]], zfhat$est, zfs$curv, 3), cbind(yfhat$x.grid[[1]], yfhat$est, yfs$curv, 5), cbind(xfhat$x.grid[[1]], xfhat$est, xfs$curv, 7))
      fs = data.frame(x=fs[,1], f=fs[,2], curv=fs[,3], lab=fs[,4])
      fs.list = rbind(fs.list, cbind(fs, which.mixt=j))
    }
  }  #end of pm loop

#---------------------------------------
  
  dens = dens.list
  dens.cont = dens.cont.list
  if (add.feature) fs = fs.list

  levplot = xyplot(y ~ x | lab, data=dens, region=FALSE, as.table=TRUE, 
  horizontal=TRUE,subscripts=TRUE, layout=c(3,3),
     scales=list(x=list(alternating=c(2,3,1), tck=0), 
     y=list(alternating=c(1,3,2), tck=0)), xlab="", ylab="",
     panel=function(x, y, subscripts, ...) {

       addSignifFeatureRegion.1d.panel = function(featureCol)
       {
         dest = list(x.grid=list(fs2$x), est=fs2$f)
         SignifFeature = as.logical(fs2$curv)
         SignifFeature.inds = which(SignifFeature)
         SGlen = length(SignifFeature.inds)
         diff.vec = diff(SignifFeature.inds)
         jump.inds = (1:length(diff.vec))[diff.vec!=1]
         num.jumps = length(jump.inds)
         
         if (num.jumps==0) 
	 panel.lines(dest$x.grid[[1]][SignifFeature.inds], dest$est[SignifFeature.inds], col=featureCol, lwd=lwd)
         
         if (num.jumps>0)
         {
           curr.inds = SignifFeature.inds[1:jump.inds[1]]
           panel.lines(dest$x.grid[[1]][curr.inds],dest$est[curr.inds], col=featureCol, lwd=lwd)
           if (num.jumps>1) 
           { 
             for (j in 2:length(jump.inds))
             {
               curr.inds = SignifFeature.inds[(jump.inds[j-1]+1):jump.inds[1]]
               panel.lines(dest$x.grid[[1]][curr.inds],dest$est[curr.inds], col=featureCol, lwd=lwd)
             }
           }
           curr.inds = SignifFeature.inds[(max(jump.inds)+1):SGlen]
           panel.lines(dest$x.grid[[1]][curr.inds],dest$est[curr.inds], col=featureCol, lwd=lwd)
         }
       }  #end of add???() function
       
       for (j in 1:length(params))
        {
	  subscripts.temp = subscripts[subscripts %in% which(dens$which.mixt==j)]
          num = as.numeric(na.omit(unique(dens$lab[subscripts.temp]))) 
          
          if (num==3 | num==5 | num==7)
          {
            xrange = current.panel.limits()$xlim
            yrange = current.panel.limits()$ylim

            xrange[1] = xrange[1]+0.1*abs(xrange[1])
            xrange[2] = xrange[2]-0.1*abs(xrange[2])

            yrange[1] = yrange[1]+0.1*abs(yrange[1])
            yrange[2] = yrange[2]-0.1*abs(yrange[2])

            densf = dens$f[subscripts.temp]

            if (add.feature) fs2 = fs[fs$which.mixt==j & fs$lab==num,]
            
            if (marginal.scale==0)
	    {
              axis.labels.orig = pretty(dens.list$f)
              axis.labels = (axis.labels.orig-min(dens.list$f))/(max(dens.list$f)-min(dens.list$f)) 
              densf = (densf-min(dens.list$f))/(max(dens.list$f)-min(dens.list$f))
              if (add.feature) fs2$f = (fs2$f-min(dens.list$f))/(max(dens.list$f)-min(dens.list$f))
            }
	    else
	    {
	       axis.labels.orig = pretty(dens.list$f[subscripts], n=4)
	       axis.labels = (axis.labels.orig-min(dens.list$f[subscripts]))/(max(dens.list$f[subscripts])-min(dens.list$f[subscripts]))
 	       densf = (densf-min(dens.list$f[subscripts]))/(max(dens.list$f[subscripts])-min(dens.list$f[subscripts]))
               if (add.feature) fs2$f = (fs2$f-min(dens.list$f[subscripts]))/(max(dens.list$f[subscripts])-min(dens.list$f[subscripts]))
	    }   

            densf = yrange[2]*densf + (1-densf)*yrange[1]

            if (add.feature) fs2$f = yrange[2]*fs2$f + (1-fs2$f)*yrange[1]

            axis.labels = yrange[2]*axis.labels + (1-axis.labels)*yrange[1]
            
            panel.lines(x=dens$x[subscripts.temp],y=densf, col=marginal.col[j], lwd=lwd)

            panel.text(x=mean(xrange), y=yrange[2]+0.05*abs(yrange[2]), labels=rev(var.label)[(num-1)/2], col=marginal.ann.col)

	    if (add.feature) addSignifFeatureRegion.1d.panel(featureCol=feature.col[j])

	    if (j==1)
 	    {
	       if (num==3) panel.axis(side="right", at=axis.labels, labels=axis.labels.orig, draw.labels=TRUE, half=FALSE, check.overlap=TRUE, line.col=marginal.ann.col, text.col=marginal.ann.col)
               if (num==7 | num==5) panel.axis(side="left", at=axis.labels, labels=axis.labels.orig, draw.labels=TRUE, half=FALSE, check.overlap=TRUE, line.col=marginal.ann.col, text.col=marginal.ann.col)
	    }	  
	  } #end if num
          else
          {
contlev = c(dens.cont[dens.cont$which.mixt==j & dens.cont$lab==num,1], 1.01*max(dens$f))

panel.contourplot(x=dens$x,y=dens$y, z=dens$f, subscripts=subscripts.temp, contour=TRUE, region=FALSE, at=contlev, col=line.col[j], lwd=lwd)

if (num==1) {
panel.axis(side="left", draw.labels=FALSE, half=FALSE); 
panel.axis(side="top", draw.labels=FALSE, half=FALSE)}
            
if (num==2) panel.axis(side="top", draw.labels=FALSE, half=FALSE)
            if (num==4) panel.axis(side="left", draw.labels=FALSE, half=FALSE)
            if (num==6) panel.axis(side="right", draw.labels=FALSE, half=FALSE)
            if (num==8) panel.axis(side="bottom", draw.labels=FALSE, half=FALSE)
if (num==9) {
panel.axis(side="right", draw.labels=FALSE, half=FALSE); 
panel.axis(side="bottom", draw.labels=FALSE, half=FALSE)}               
          } #end elseif
        } #end of j loop
      } # end of panel function
      , ...)
                  
    plot(levplot)
}


