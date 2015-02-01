library("MASS")
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
GPmatcomp<-function(p,n,r,SigmaA,SigmaE,k,directory,cor, graph) {
  
  LisG<-vector("list",k)
  LisP<-vector("list",k)
  varmat<-matrix(0,nrow=k,ncol=p)
  Peig<-matrix(0,nrow=k,ncol=p)
  Gvarmat<-matrix(0,nrow=k,ncol=p)
  newmx<-matrix(0,nrow=k,ncol=p)
  remssl<-matrix(0,nrow=k,ncol=p)
  for (i in 1:k){
    Y <- gendat(p,n,r,SigmaA,SigmaE)
    Yssp <- ssp(Y,n,r)
    remssl[i,]<-reml(Yssp$W,Yssp$B,n,r)$lamB
    LisP[[i]]<-var(Y)
    names(LisP)[i]<-paste(i,".csv",sep="")
    if (cor==TRUE) {
      L<-Y
      L<-scale(L, center=TRUE, scale=TRUE)
      LisP[[i]]<-cor(L)
    }
    LisG[[i]]<-manov(Yssp$W, Yssp$B, n,r)$G
    names(LisG)[i]<-paste(i,".csv",sep="")
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

  if (graph==TRUE) {
  
  mypath <- file.path(paste(directory),paste("PGP","_boxplot" ,paste("_p",p,"n",n,"r",r,"k",k,sep=""), ".pdf", sep = ""))
  pdf(file=mypath)
  boxplot(ProjPG,col="grey",boxwex=0.4,xlab="Eigenvector",ylab="Projection of P through G")
  abline(0,0,lty=2)
  dev.off()
  
  ###boxplot for eigenvaluees and variance of P matrix
  mypath <- file.path(paste(directory),paste("Peig&Var","_boxplot" ,paste("_p",p,"n",n,"r",r,"k",k,sep=""), ".pdf", sep = ""))
  pdf(file=mypath)
  boxplot(Peig,col="grey",boxwex=0.5,xaxis=NULL,xlab="Trait or Eigenvector",ylab="Phenotypic variance")
  abline(mean(diag(SigmaE)+diag(SigmaA)),0,lty=2)
  dev.off()
  
  ###boxplot for eigenvaluees and variance of G matrix
  mypath <- file.path(paste(directory),paste("Geig&Var","_boxplot" ,paste("_p",p,"n",n,"r",r,"k",k,sep=""), ".pdf", sep = ""))
  pdf(file=mypath)
  boxplot(newmx,col="grey",boxwex=0.5,xaxis=NULL,xlab="Trait or Eigenvector",ylab="Genetic variance")
  abline(mean(diag(SigmaA)),0,lty=2)
  dev.off()
  
  mypath <- file.path(paste(directory),paste("Peig&Geig","_boxplot" ,paste("_p",p,"n",n,"r",r,"k",k,sep=""), ".pdf", sep = ""))
  pdf(file=mypath)
  par(mfrow=c(1,2))
  boxplot(Peig,col="grey",boxwex=0.5,xaxis=NULL,xlab="Trait or Eigenvector",ylab="Phenotypic variance")
  abline(mean(diag(SigmaE)+diag(SigmaA)),0,lty=2)
  boxplot(newmx,col="grey",boxwex=0.5,xaxis=NULL,xlab="Trait or Eigenvector",ylab="Genetic variance")
  abline(mean(diag(SigmaA)),0,lty=2)
  dev.off()

  #### implements spectral distribution as given in Iain's Code
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
  }
  else if (missing(graph)){
    
  }
  return(list(LisP=LisP,LisG=LisG,remssl=remssl,Peig=Peig,newmx=newmx))
}


TWshift<-function(vec,n,p) { ### n = no individuals, p = no traits
  m<-rep(0,length(vec))
  a<-((sqrt(p)+sqrt(n))^(2))
  b<-((sqrt(n)+sqrt(p))*((1/sqrt(p)+1/sqrt(n))^(1/3)))
  for (i in 1:length(vec)) {
    m[i]<-(n*vec[i]-a)/b
  }
  return (m)
}

if (FALSE) {
n<-sqrt(q$newmx[,1])/mean(q$newmx[,1])
l<-sqrt(q$Peig)/colMeans(q$Peig)
l2<-sqrt(q$newmx[,1:3])/rep(colMeans(q$newmx[,1:3]),each=1000)
n1<-sqrt(q$newmx[,1])/mean(q$newmx[,1])
n2<-sqrt(q$newmx[,2])/mean(q$newmx[,2])
n3<-sqrt(q$newmx[,3])/mean(q$newmx[,3])

}

if (FALSE) {
z<-data.frame(matrix(NA,nrow=1000,ncol=3))
z[,1]<-n1
z[,2]<-n2
z[,3]<-n3
mypath <- file.path(paste(directory),paste("CV" ,paste("_p",10,"n",50,"r",5,"k",1000,sep=""), ".pdf", sep = ""))
pdf(file=mypath)
boxplot(l2,boxwex=0.3, col='grey', xaxis=NULL, ylim=c(0.8,max(l2)),xlab="Eigenvector",ylab="coefficient of variance")
boxplot(l[,1:3],boxwex=0.3, col='white', xaxis=NULL,add=TRUE,at=1:3+0.3)
legend(x="bottomleft", fill=c("grey","white"), legend=c("CV(G)","CV(P)"))
dev.off()

rand<-GPmatcomp(10,50,5,0*diag(10),diag(10),1000, "c:/ABakeSumProj",cor=TRUE)
ranTW<-TWshift(rand$Peig[,1],250,10)
ransimTW<-rtw(1000,beta=1)
ransimTW2<-rtw(1000,beta=2)

mypath <- file.path(paste("C:/ABakeSumProj"),paste("TW_sim_versus_obs2" ,paste("_p",10,"n",50,"r",5,"k",1000,sep=""), ".pdf", sep = ""))
pdf(file=mypath)
boxplot(ransimTW,boxwex=0.2, col='grey', xaxis=NULL,at=1)
boxplot(ranTW,boxwex=0.2, col='white', xaxis=NULL,add=TRUE,at=0.8)
boxplot(ransimTW,boxwex=0.2, col='grey42', xaxis=NULL,add=TRUE,at=1.2)
legend(x="bottomleft", fill=c("white","grey","grey42"), legend=c("shifted","simulatedTW1","simulatedTW2"))
dev.off()

mypath <- file.path(paste("C:/ABakeSumProj"),paste("QQTW" ,paste("_p",10,"n",50,"r",5,"k",10000,sep=""), ".pdf", sep = ""))
pdf(file=mypath)
qqplot(ransimTW,ranTW,xlab="Shifted and rescaled Correlation matrix eigenvalues",ylab="Randomly simulated Tracy-Widom rescaling")
qqline(ransimTW,distribution=rtw)
dev.off()

mypath <- file.path(paste("C:/ABakeSumProj"),paste("Coefficient_of_variation" ,paste("_p",10,"n",50,"r",5,"k",10000,sep=""), ".pdf", sep = ""))
pdf(file=mypath)
plot(1:10,sqrt(diag(var(Nrand$newmx)))/colMeans(Nrand$newmx),ylab="CV",xlab="Eigenvalues",pch=0,bg="grey")
points(1:10,sqrt(diag(var(Nrand$Peig)))/colMeans(Nrand$Peig),pch=15)
legend(x="topright",fill=c("white","black"), legend=c("G_eigenvalues","P_eigenvalues"))
dev.off()

mypath <- file.path(paste("C:/ABakeSumProj"),paste("Coefficient_of_variation for P" ,paste("_p",10,"n",50,"r",5,"k",10000,sep=""), ".pdf", sep = ""))
pdf(file=mypath)
plot(1:10,sqrt(diag(var(Nrand$Peig)))/colMeans(Nrand$Peig),ylab="CV",xlab="Eigenvalues",pch=15,bg="grey")
dev.off()
}

Multrun<-function(p,n,r) {
  z<-vector("list",length(p))
for ( i in 1:length(p)) {
  z[[i]]<-GPmatcomp(p[i],n[i],r[i],0*diag(p[i]),diag(p[i]),1000, "C:/ABakeSumProj", cor=TRUE)
  
  }
  return (z)
}

GreaterthanNindx<-function(index,n,dir) {
  m<-read.csv(index,header=TRUE,stringsAsFactors=FALSE)
  q<-m[m[,9]>=5,]
  write.csv(q,file=paste(dir,"nsubmats.csv",sep="/"),row.names=FALSE)
  return(q)
}

PonNconv<-function(p,r,index) {
  z<-read.csv(file=index, header=TRUE,stringsAsFactors=FALSE)
  m<-rep(0,length(z[,1]))
  for (i in 1:length(z[,1])) {
    m[i]<-z[i,17] -z[i,17]%%r
    m[i]<-m[i]/r
  }  
  return (m)
}

TWidomTest<-function(index, r,p, k,directory, graph,PaGlist,mean) { ### Issue with Pmax and TWd
  q<-vector("list",length(PaGlist[,1]))
  Shift<-vector("list",length(PaGlist[,1]))
  Q<-PonNconv(p,r,index)
  ObsTW<-matrix(0,nrow=length(PaGlist[,1]),ncol=length(PaGlist[1,]))
  Md<-matrix(0,nrow=length(PaGlist[,1]),ncol=length(PaGlist[1,]))
  Sd<-matrix(0,nrow=length(PaGlist[,1]),ncol=length(PaGlist[1,]))
  TSd<-matrix(0,nrow=length(PaGlist[,1]),ncol=length(PaGlist[1,]))
  for (i in 1:length(PaGlist[,1])) {
    q[[i]]<-GPmatcomp(p,Q[i],r,0*diag(p),diag(p),k,directory,TRUE,graph)
    Shiftsto<-matrix(0,nrow=length(q[[i]]$Peig[,1]),ncol=p)
    l<-rep(0,k)
    Pmax<-0
    TWd<-0
    for (j in 1:length(PaGlist[1,])) {
    l<-q[[i]]$Peig[,j]
    qden<-density(l)
    if (mean==FALSE) {
    Md[i,j]<-sum(((qden$x*qden$y))*(qden$x[2]-qden$x[1])) ### first raw moment = E(X)=mean
    }
    else {
      Md[i,j]<-mean(l)
    }
    TSd[i,j]<-sqrt(sum(((qden$x-Md[i,j])^2*qden$y))*(qden$x[2]-qden$x[1])) ### sqrt of the 2nd central moment about the mean
    Pmax<-PaGlist[i,j]
    TWd<-(-1.206+(1.268/TSd[i,j])*(Pmax-Md[i,j])) #### used MTW and STW^2 from Saccenti et al
    Shiftsto[,j]<-(-1.206+(1.268/TSd[i,j])*(q[[i]]$Peig[,j]-Md[i,j]))
    ObsTW[i,j]<-TWd
    }
    Shift[[i]]<-Shiftsto
  }
  return (list(ObsTW=ObsTW,Peig=l,Md=Md,TSd=TSd,Shift=Shift))
}
### Used command - TWidomTest("G:/GIThub/Rproject/nsubmats.csv",5,5,10000, "G:/GIThub/Rproject",FALSE, test2[[2]])
##### Generating Tracy-Widom plot

if (FALSE) {
##### Code for writing final version of data into folders ### Note median p/n ~ 0.012
permsto<-GPmatcomp(5,83,5,0*diag(5),diag(5), 1000,"C:/Users/s4284361/Documents/GitHub/Rproject/Simulatedstorage/p=5/Graphs", FALSE,TRUE)
WriteMatList(permsto[[1]],"C:/Users/s4284361/Documents/GitHub/Rproject/Simulatedstorage/p=5/Pmat")
WriteMatList(permsto[[2]],"C:/Users/s4284361/Documents/GitHub/Rproject/Simulatedstorage/p=5/Gmat")
setwd("C:/Users/s4284361/Documents/GitHub/Rproject/Simulatedstorage/p=5")
write.csv(permsto[[4]],row.names=FALSE,file="SimulatedPeigenvalues.csv")
write.csv(permsto[[5]],row.names=FALSE,file="SimulatedGeigenvalues.csv")

permsto<-GPmatcomp(10,166,5,0*diag(10),diag(10), 1000,"C:/Users/s4284361/Documents/GitHub/Rproject/Simulatedstorage/p=10/Graphs", FALSE,TRUE)
WriteMatList(permsto[[1]],"C:/Users/s4284361/Documents/GitHub/Rproject/Simulatedstorage/p=10/Pmat")
WriteMatList(permsto[[2]],"C:/Users/s4284361/Documents/GitHub/Rproject/Simulatedstorage/p=10/Gmat")
setwd("C:/Users/s4284361/Documents/GitHub/Rproject/Simulatedstorage/p=10")
write.csv(permsto[[4]],row.names=FALSE,file="SimulatedPeigenvalues.csv")
write.csv(permsto[[5]],row.names=FALSE,file="SimulatedGeigenvalues.csv")

}
#####
####  rm(list = setdiff(ls(), lsf.str())) removes all non-function arguments from environment


if (FALSE) {
  load("F:/randomised_empirical_shift.RData")
  library(RMTstat)
  RTW<-rtw(10000,1)
  ### Graph comparing observed data to theoretical distribution.   
  m<-hist(RTW,breaks=50, xlab="TW statistic",main=NULL,freq=FALSE,ylim=c(0,0.38))
  abline(v=qtw(0.05,lower.tail=FALSE),lty=2)
  mx<-max(m$density)
  mn<-min(m$density)
  sub<-mx-mn
  for (i in 1:length(teststo[[1]])) {
    points(type="l",x=c(teststo[[1]][i],teststo[[1]][i]),y=c(mx+0.05*sub,mx+0.15*sub))
  }
  
  
  ### Graph comparing fit of null data to theoretical distribution
  masto<-matrix(0,nrow=430000,ncol=5)
  for (i in 1:length(teststo$Shift)) {
    masto[((10000*(i-1)+1):(10000*i)),1:5]<-teststo$Shift[[i]][,1:5]
  }
  
  mypath=
  pdf(file=mypath)
  hist(RTW,breaks=50,freq=FALSE,ylim=c(0,0.35), main=NULL,xlab="Random Tracy-Widom values")
  lines(density(masto[,1]),lwd=2,col="black")
  lines(density(masto[,3]),lwd=2,col="#E69F00")
  lines(density(masto[,5]),lwd=2, col="#0072B2")
  legend(x="topright", lty=c(1,1,1),col=c("black","#E69F00","#0072B2"), legend=c(expression(paste("Shifted ",lambda[1],sep="")),expression(paste("Shifted ",lambda[3],sep="")),expression(paste("Shifted ",lambda[5],sep=""))),bty="n")
  dev.off()
  
}


