
do.features<- function(retpath='.',p=4,g=2,distr="mst",
clusters=c("base","mould"),
data="0min",panel="panel1",group="LNP_hi",
samples,markers=c("p-SFK","BCL2" ,"CD20","p-ERK"),class="jcm") {

#sample names are divided into two groups according to the "info.csv" file.
#or just colbind their feature matrices.

# distribution parameters
#parameters <- c("error","loglik","bic","pro","mu","sigma",
#"dof","delta","mode","scale","shape","orientation")

distr<- tolower(distr)

if(is.na(match(distr,c("mvn","mvt","msn","mst"))))
stop("distribution specified is not available!")

# for symmetric matrix sigma

crossnames1 <- rep('',p*(p+1)/2)
counter=1
for(j in 1:p) {
for(i in j:p) {
crossnames1[counter] <- paste(markers[i],markers[j],sep='.')
counter = counter + 1
}}

# for orientation

crossnames2 <- rep('',p*p)
counter=1
for(j in 1:p) {
for(i in 1:p) {
crossnames2[counter] <- paste(markers[i],markers[j],sep='.')
counter = counter + 1
}}

#----------------------------------------------

rnames <- c("error","loglik","bic")

for(h in 1:g) {
rnames <- c(rnames,paste("pro",clusters[h],sep='.'))
}

for(h in 1:g) {
rnames <- c(rnames,
paste(paste("mu",clusters[h],sep='.'),markers,sep='.'))
}
for(h in 1:g) {
rnames <- c(rnames,
paste(paste("sigma",clusters[h],sep='.'),crossnames1,sep='.'))}

for(h in 1:g) {
rnames <- c(rnames,paste("dof",clusters[h],sep='.'))
}

for(h in 1:g) {
rnames <- c(rnames,
paste(paste("delta",clusters[h],sep='.'),markers,sep='.'))}

for(h in 1:g) {
rnames <- c(rnames,
paste(paste("mode",clusters[h],sep='.'),markers,sep='.'))}

for(h in 1:g) {
rnames <- c(rnames,
paste("scale",clusters[h],sep='.'))}

for(h in 1:g) {
rnames <- c(rnames,
paste(paste("shape",clusters[h],sep='.'),markers,sep='.'))}

for(h in 1:g) {
rnames <- c(rnames,
paste(paste("orientation",clusters[h],sep='.'),crossnames2,sep='.'))
}

fmat <- array(0,c(length(rnames),length(samples)))


for(file in samples) {

if(tolower(class)=="jcm")
ret <- dget(paste(retpath,'/',data,'_',panel,'_',file,'_',distr,'_g=',g,'_jcm.ret',sep=''))
else
ret <- dget(paste(retpath,'/',data,'_',panel,'_',file,'_',distr,'_g=',g,'.ret',sep=''))


if(tolower(class)=="jcm") {

#the if statement is for compatibility to old versions

if(is.null(ret$modpts)) {

ret$sigma2 <- doentro(p,g, distr, class="jcm", ret) 
ret$modpts <- emmixMOD(p, g, distr, ret$mu, ret$sigma2, 
ret$dof, ret$delta, nrand = 10000)
}

}

#----------------------

#specific to BCR data

if(!is.na(match("base",tolower(clusters))))
{

obj <- ret
# mould
if( sum( (obj$modpts[2:3,1])^2 ) > sum( (obj$modpts[2:3,2])^2 ))
ret <- doswap(obj,c(2,1),class)
obj <- NULL
}

#---------------------


# add some new charatertics

obj <- sigma2decomp(ret$sigma)[[2]]
ret$scale <- obj$scale
ret$shape <- obj$shape
ret$orientation <- obj$orientation


fmat[,match(file,samples)] <- c(ret$error,ret$loglik,
ret$bic,ret$pro,ret$mu,
lowertriangle(ret$sigma,p,g),ret$dof,ret$delta,
ret$modpts,ret$scale,ret$shape,ret$orientation)
}

dimnames(fmat) <- list(rnames,samples)

#write to gct format
gct <- list()
gct$data <- fmat

if(tolower(class)=="jcm") {
filename1 <- paste('batch','features',data,panel,
group,distr,paste("g=",g,sep=''),"jcm.txt",sep='_')
filename2 <- paste('batch','features',data,panel,
group,distr,paste("g=",g,sep=''),"jcm.gct",sep='_')
}
else {
filename1 <- paste('batch','features',data,panel,
group,distr,paste("g=",g,".txt",sep=''),sep='_')
filename2 <- paste('batch','features',data,panel,
group,distr,paste("g=",g,".gct",sep=''),sep='_')
}


write.table(fmat, filename1,sep='\t')
write.gct(gct, filename2)

}

#-----------------------------------------



write.gct <-
#
# save a GCT result to a file, ensuring the filename has the extension .gct
#
function(gct, filename)#, check.file.extension=TRUE)
{


# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2003-2006) by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.

# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for its
# use, misuse, or functionality.




#if(check.file.extension) {
#		filename <- check.extension(filename, ".gct") 
#	}

#	if(is.null(gct$data)) {
#		exit("No data given.")
#	}
#	if(is.null(row.names(gct$data))) {
#		exit("No row names given.")
#	}
#	if(is.null(colnames(gct$data))) {
#		exit("No column names given.")
#	}
	
	rows <- dim(gct$data)[1]
	columns <- dim(gct$data)[2]
	
	if(rows!=length(row.names(gct$data))) {
		stop("Number of data rows (", rows, ") not equal to number of row names (", length(row.names(gct$data)), ").")
	}
	if(columns!=length(colnames(gct$data))) {
		stop("Number of data columns (", columns , " not equal to number of column names (", length(colnames(gct$data)), ").")
	}
		
	if(!is.null(gct$row.descriptions)) {
		if(length(gct$row.descriptions)!=rows) {
			stop("Number of row descriptions (", length(gct$row.descriptions), ") not equal to number of row names (", rows, ").")
		}
	}
	
	row.descriptions <- gct$row.descriptions
	if(is.null(row.descriptions)) {
	    row.descriptions <- ''
	}
	
	m <- cbind(row.names(gct$data), row.descriptions, gct$data)
	
	f <- file(filename, "w")
	
	on.exit(close(f))
	
	cat("#1.2", "\n", file=f, append=TRUE, sep="")
	cat(rows, "\t", columns, "\n", file=f, append=TRUE, sep="")
	cat("Name", "\t", file=f, append=TRUE, sep="")
	cat("Description", file=f, append=TRUE, sep="")
	names <- colnames(gct$data)
	
	for(j in 1:length(names)) {
		cat("\t", names[j], file=f, append=TRUE, sep="")
	}

	cat("\n", file=f, append=TRUE, sep="")
	write.table(m, file=f, append=TRUE, quote=FALSE, sep="\t", eol="\n", col.names=FALSE, row.names=FALSE)
	return(filename)
}



sigma2decomp<-function (sigma, G = NULL, tol = NULL, ...) 
{

#Warning! 
#The function sigma2decomp() is a copy from mclust;
#please evoke mclust and read the copyright issue.

    
    dimSigma <- dim(sigma)
    if (is.null(dimSigma)) 
        stop("sigma improperly specified")
    d <- dimSigma[1]
    if (dimSigma[2] != d) 
        stop("sigma improperly specified")
    l <- length(dimSigma)
    if (l < 2 || l > 3) 
        stop("sigma improperly specified")
    if (is.null(G)) {
        if (l == 2) {
            G <- 1
            sigma <- array(sigma, c(dimSigma, 1))
        }
        else {
            G <- dimSigma[3]
        }
    }
    else {
        if (l == 3 && G != dimSigma[3]) 
            stop("sigma and G are incompatible")
        if (l == 2 && G != 1) 
            sigma <- array(sigma, c(d, d, G))
    }
    decomp <- list(d = d, G = G, scale = rep(0, G), shape = matrix(0, 
        d, G), orientation = array(0, c(d, d, G)))
    for (k in 1:G) {
        ev <- eigen(sigma[, , k], symmetric = TRUE)
        temp <- log(ev$values)
        logScale <- sum(temp)/d
        decomp$scale[k] <- exp(logScale)
        decomp$shape[, k] <- exp(temp - logScale)
        decomp$orientation[, , k] <- ev$vectors
    }
    if (is.null(tol)) 
        tol <- sqrt(.Machine$double.eps)
    scaleName <- "V"
    shapeName <- "V"
    orientName <- "V"
    uniq <- function(x, tol = sqrt(.Machine$double.eps)) {
        abs(max(x) - min(x)) < tol
    }
    if (uniq(decomp$scale)) {
        decomp$scale <- decomp$scale[1]
        scaleName <- "E"
    }
    if (all(apply(decomp$shape, 1, uniq, tol = tol))) {
        decomp$shape <- decomp$shape[, 1]
        if (all(uniq(decomp$shape, tol = tol))) {
            shapeName <- "I"
            decomp$shape <- rep(1, d)
        }
        else {
            shapeName <- "E"
        }
    }
    if (all(apply(matrix(decomp$orientation, nrow = d * d, ncol = G), 
        1, uniq, tol = tol))) {
        decomp$orientation = decomp$orientation[, , 1]
        if (all(apply(cbind(decomp$orientation, diag(d)), 1, 
            uniq, tol = tol))) {
            orientName <- "I"
            decomp$orientation <- NULL
        }
        else {
            orientName <- "E"
        }
    }
    modelName <- paste(c(scaleName, shapeName, orientName), collapse = "")
    c(list(modelName = modelName, decomp))
}




lowertriangle <- function(a,p,g) {

a <- array(a,c(p,p,g))

element <- rep(0,g*(p+1)*p/2)

counter=1

for(h in 1:g) {

for(j in 1:p) {

for(i in j:p) {

element[counter] <- a[i,j,h]

counter = counter + 1

}
}
}

element
}



doswap <- function(obj,st=c(2,1),class='') {

g <- length(obj$pro)

if(length(st)!=g) stop("length of st does not match!")

ret <- obj

for(h in 1:g)
{
ret$pro[h]    <- obj$pro[st[h]]
ret$mu[,h]    <- obj$mu[,st[h]]
ret$sigma[,,h]<- obj$sigma[,,st[h]]
ret$dof[h]    <- obj$dof[st[h]]
ret$delta[,h] <- obj$delta[,st[h]]
ret$modpts[,h]<- obj$modpts[,st[h]]

if(tolower(class)=="jcm")
{
ret$theta[,h] <- obj$theta[,st[h]]
ret$thetu[h]  <- obj$thetu[st[h]]
}

}

ret
}

