
library("MASS")
library("Matrix")
###Start of Iain's code#####

genspec <- function(c1,c2,x,range,step) {
  # compute limiting density of eigenvalues of H - xE for n_1 H ~
  # W_p(n_1,I) independent of n_2 E ~ W_p(n_2,I) with c1 = lim p/n,
  # c2 = lim p/n2.  range = c(lo,hi) is range of z-values on real
  # axis with stepwidth step
  
  b1 <- 1+c2*x - (x+c1)
  b2 <- -x*(c1 + c2 - c1*c2)
  
  argv <- seq( range[1], range[2], by = step)
  npts <- length(argv)
  dens <- rep(0,npts)
  
  for (i in  1:npts) {
    z <- argv[i]
    a0 <- -1
    a1 <- -z + b1
    a2 <- -c1*z + c2*x*z + b2
    a3 <- c1*c2*z*x
    av <- c(a0,a1,a2,a3)
    dens[i] <- max(Im(polyroot(av)))/pi
  }
  list( argv=argv, dens=dens)
}


manov <- function(W,B,n,r){
  # compute reml estimates Bhat and What of effect and error
  # covariance matrices and eigenvalues lamB of Bhat. Balanced case,
  # n = no. groups, r = no. reps per group
  # method avoids evals of W^{-1}B
  
  return(list(eige1=(eigen(B-W)$values/r),G=(B-W)/r))
}

ssp <- function(Y,n,r){
  # compute within-group and between-group mean square matrices for
  # balanced case, n = no. groups, r = no. reps per group
  
  gp <- rep(1:n, each=r)
  p <- dim(Y)[2]
  
  Ymeans <- matrix(0, nrow = n, ncol = p)
  E      <- matrix(0, nrow = p, ncol = p)
  
  for (i in 1:n){
    Yt <- Y[gp==i,]
    Ymeans[i,] <- colMeans(Yt)
    E <- E + (r-1)*cov(Yt)
  }
  nE <- n*(r-1)
  return( list(W=E/nE, B = r*cov(Ymeans)))
}

reml <- function(W,B,n,r){
  # compute reml estimates Bhat and What of effect and error
  # covariance matrices and eigenvalues lamB of Bhat. Balanced case,
  # n = no. groups, r = no. reps per group
  # method avoids evals of W^{-1}B
  
  L <- solve(chol(W))
  
  gen <- eigen( t(L) %*% B %*% L, symmetric= TRUE)
  lamhat <- gen$values
  P <- solve(t(L)) %*% gen$vectors
  
  tlam <- ifelse( lamhat > 1, lamhat-1, 0)
  Bhat <- P %*% diag(tlam) %*% t(P) /r
  lamB <- eigen( Bhat, only.values=TRUE)$values
  
  nE <- n*(r-1)
  nH <- n-1
  What <- ( nH*(B - r*Bhat) + nE*W )/(nH + nE)
  list( Bhat=Bhat, What=What, lamB=lamB)
}


## end of Iain's code

##New generate data funtion which can handle case with no random effect for line
## modified from code supplied by Iain (personal communication)
gendat <- function(p,n,r, SigmaA, SigmaE){
  # generate balanced one-way MANOVA with random effects covariance
  # SigmaA and error covariance SigmaE
  # p-dimensional observations on n groups with r replicates per group
  
  gp <- rep(1:n, each=r)
  zero <- rep(0, p)
  if (all(diag(SigmaA)==zero)) {
    error <- mvrnorm(n*r, zero, SigmaE)
    return(gendat=error)
  } 
  else {
    randef <- mvrnorm(n, zero, SigmaA)
    vcpt <- randef[gp,]       # repeat each random effect r times
    error <- mvrnorm(n*r, zero, SigmaE)
    return( gendat = vcpt + error)
  }
  
}
###### Function which can take p=no. traits, n=no. lines, r = no. reps in each line, SigmaA=line randomeffect covariance matrix, SigmaE=error covariance matrix,k=no. simulated datasets
###### simulated data is made a multivariate normal with covariance matrix equal to either SigmaA or SigmaE to get differ
###### Directory is the location wher you wish the graphs to be saved. It needs ti be done as a string. Eg if you wish to save to C:/testfolder you must input "C:/testfolder"
## Start of work by Alex
GPmatcomp<-function(p,n,r,SigmaA,SigmaE,k,directory) {
  
  LisG<-vector("list",k)
  LisP<-vector("list",k)
  LisGR<-vector("list",k)
  varmat<-matrix(0,nrow=k,ncol=p)
  Peig<-matrix(0,nrow=k,ncol=p)
  Gvarmat<-matrix(0,nrow=k,ncol=p)
  newmx<-matrix(0,nrow=k,ncol=p)
  remssl<-matrix(0,nrow=k,ncol=p)
  NePD<-matrix(0,nrow=k,ncol=p)
  for (i in 1:k){
    Y <- gendat(p,n,r,SigmaA,SigmaE)
    Yssp <- ssp(Y,n,r)
    remssl[i,]<-reml(Yssp$W,Yssp$B,n,r)$lamB
    LisP[[i]]<-var(Y)
    LisP[[i]]<-cov2cor(LisP[[i]])
    LisG[[i]]<-manov(Yssp$W, Yssp$B, n,r)$G
    NePD[i,]<-eigen(nearPD(LisG[[i]],corr=TRUE)$mat)$values
    LisGR[[i]]<-nearPD(LisG[[i]],corr=TRUE)$mat
    newmx[i,] <- manov(Yssp$W, Yssp$B, n,r)$eige1
    Gvarmat[i,]<-sort(diag(manov(Yssp$W, Yssp$B, n,r)$G),decreasing=TRUE)
    varmat[i,]<-sort(diag(var(Y)),decreasing=TRUE)
    Peig[i,]<-eigen(LisP[[i]])$values
  }
  
  
  ProjPG<-matrix(0,nrow=k,ncol=p)
  for (i in 1:k) {
    for (j in 1:p){
      ProjPG[i,j]<-t((eigen(LisP[[i]],symmetric=TRUE)$vectors[,j]))%*%LisG[[i]]%*%(eigen(LisP[[i]],symmetric=TRUE)$vectors[,j])
    }
  }
    
  CorProjPG<-matrix(0,nrow=k,ncol=p)
  for (i in 1:k) {
    for (j in 1:p){
      CorProjPG[i,j]<-as.vector(t((eigen(LisP[[i]],symmetric=TRUE)$vectors[,j]))%*%LisGR[[i]]%*%(eigen(LisP[[i]],symmetric=TRUE)$vectors[,j]))
    }
  }
  mypath <- file.path(paste(directory),paste("CorPGP","_boxplot" ,paste("_p",p,"n",n,"r",r,"k",k,sep=""), ".pdf", sep = ""))
  pdf(file=mypath)
  boxplot(CorProjPG,col="grey",boxwex=0.4,xlab="Eigenvector",ylab="Projection of P through G")
  abline(0,0,lty=2)
  dev.off()
  
  mypath <- file.path(paste(directory),paste("PGP","_boxplot" ,paste("_p",p,"n",n,"r",r,"k",k,sep=""), ".pdf", sep = ""))
  pdf(file=mypath)
  boxplot(ProjPG,col="grey",boxwex=0.4,xlab="Eigenvector",ylab="Projection of P through G")
  abline(0,0,lty=2)
  dev.off()
  
  z<-ProjPG+1
  
  mypath <- file.path(paste(directory),paste("PGP_cor&shiftedvar","_boxplot" ,paste("_p",p,"n",n,"r",r,"k",k,sep=""), ".pdf", sep = ""))
  pdf(file=mypath)
  boxplot(CorProjPG,col="white",boxwex=0.4,xlab="Eigenvector",ylab="Projection of P through G", ylim=c(0.85,1.15))
  boxplot(z,add=TRUE, col="grey",at=1.5:(p+0.5), boxwex=0.4,xaxt='n')
  legend(x="bottomright", c("Pcor'GcorPcor","(P'GP+1)"), fill=c("white","grey"))
  abline(1,0,lty=2)
  dev.off()
  
  q<-Peig-NePD
  mypath <- file.path(paste(directory),paste("DifeigPG","_boxplot" ,paste("_p",p,"n",n,"r",r,"k",k,sep=""), ".pdf", sep = ""))
  pdf(file=mypath)
  boxplot(q,col="white",boxwex=0.4,xlab="Eigenvector",ylab="eig(P)-eig(G)")
  legend(x="topright", c("eig(P)-eig(G)"), fill=c("white"))
  abline(1,0,lty=2)
  dev.off()
    
  
  ###boxplot for eigenvaluees and variance of P matrix
  mypath <- file.path(paste(directory),paste("Peig&Var","_boxplot" ,paste("_p",p,"n",n,"r",r,"k",k,sep=""), ".pdf", sep = ""))
  pdf(file=mypath)
  boxplot(Peig,col="grey",boxwex=0.4,xaxis=NULL,xlab="Trait or Eigenvector",ylab="Phenotypic variance")
  boxplot(varmat,add=TRUE, at=1.5:(p+0.5), boxwex=0.4,xaxt='n')
  abline(mean(diag(SigmaE)+diag(SigmaA)),0,lty=2)
  dev.off()
  
  ###boxplot for eigenvaluees and variance of G matrix
  mypath <- file.path(paste(directory),paste("Geig&Var","_boxplot" ,paste("_p",p,"n",n,"r",r,"k",k,sep=""), ".pdf", sep = ""))
  pdf(file=mypath)
  boxplot(newmx,col="grey",boxwex=0.4,xaxis=NULL,xlab="Trait or Eigenvector",ylab="Genetic variance")
  boxplot(Gvarmat,add=TRUE, at=1.5:(p+0.5), boxwex=0.4,xaxt='n')
  abline(mean(diag(SigmaA)),0,lty=2)
  dev.off()
  # Peig and varmat not calculated using REML, newmx is REML estimate of eigenvalues
  
  #new REML estimates of G var matrix
  
  # Generates eigenvalues for 200 simulated G matrices.
  #need to extract unconstrained REML estimates of variances for each trait.
  
  # Need to generate density histograms to see if marchenko-Pastur distribution fits
  
  #Values used from aug13 to generate limiting spectral distribution of eigenvalues for G
  #No. lines
  # no replicates per line
  # no. traits measured
  nH <- n-1
  nE <- n*(r-1)
  c1 <- p/nH
  c2 <- p/nE
  x <- 1
  range <- c(-2,2)
  step <- .05
  
  spec <- genspec(c1,c2,x,range,step)
  
  mypath <- file.path(paste(directory),paste("Geig","_Hist" , paste("_p",p,"n",n,"r",r,"k",k,sep=""), ".pdf", sep = ""))
  pdf(file=mypath)
  hist(newmx,breaks=30, freq=FALSE,xlab="Eigenvalues",main='') # HISTOGRAM of G matrix values
  lines(spec$argv/5,5*spec$dens)
  dev.off()
  
  mypath <- file.path(paste(directory),paste("Peig","_Hist" ,paste("_p",p,"n",n,"r",r,"k",k,sep=""), ".pdf", sep = ""))
  pdf(file=mypath)
  hist(Peig,breaks=30, freq=FALSE,xlab="Eigenvalues") #histogram of Pmatrix values
  dev.off()
  
  
  mypath <- file.path(paste(directory),paste("REMLvsUncon","_boxplot" ,paste("_p",p,"n",n,"r",r,"k",k,sep=""), ".pdf", sep = ""))
  pdf(file=mypath)
  boxplot(newmx,col="grey",boxwex=0.4,xaxis=NULL,xlab="Trait or Eigenvector",ylab="Genetic Variance")
  boxplot(remssl,add=TRUE, at=1.5:(p+0.5), boxwex=0.4,xaxt='n')
  abline(mean(diag(SigmaE)+diag(SigmaA)),0,lty=2)
  dev.off()

  mypath <- file.path(paste(directory),paste("NePDvsUncon","_boxplot" ,paste("_p",p,"n",n,"r",r,"k",k,sep=""), ".pdf", sep = ""))
  pdf(file=mypath)
  boxplot(newmx,col="grey",boxwex=0.4,xaxis=NULL,xlab="Trait or Eigenvector",ylab="Genetic Variance")
  boxplot(NePD,add=TRUE, at=1.5:(p+0.5), boxwex=0.4,xaxt='n')
  abline(mean(diag(SigmaE)+diag(SigmaA)),0,lty=2)
  dev.off()
  
  return(list(LisP=LisP,LisG=LisG,remssl=remssl,LisGR=LisGR,NePD=NePD,Peig=Peig))
}




