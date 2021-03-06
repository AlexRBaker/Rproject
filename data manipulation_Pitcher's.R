
if (FALSE) { ## Wrapping earlier code in if statement so it doesn't get run by accident.
######### Creating suitable data frames for storinf the data

p<-list.files(path="G:/GIThub/Pitchers_PTRS2014/Data/Gmats_&_means_as_CSVs")
q<-list.files(path="G:/GIThub/Pitchers_PTRS2014/Data/Gmats_Cor_as_CSVs")
z<-read.csv("G:/GIThub/Pitchers_PTRS2014/Data/MatrixIndexFinal.csv",header=TRUE)
datdes<-data.frame(matrix(NA,nrow=(length(p)+length(q)),ncol=15))
colnames(datdes)<-c("Paper","Organism","Raw data avilable", "Gmatrix","Pmatrix","cormatrix","directory","data link","Reference","Taxon1","Taxon2")

##Indicator function which creates no. row entries with a 0 or 1 for the Genetic covariance/correlation matrix. 1=presence, 0=absence with a corresponding paper key
for (i in 1:(length(p)+length(q))) {
  if (i<=length(p)) {
    datdes[i,1]<-p[i]
    datdes[i,4]<-1
    datdes[i,6]<-0
  }  else {
   datdes[i,1]<-q[i-length(p)]
   datdes[i,4]<-0
   datdes[i,6]<-1
  }
  
}

dades2<-datdes[!grepl("list", datdes[,1]),]   ### subset of data with only csv titles from directory files

substr(dades2[1,1],1,7) ## Checking function of substring command

##### Retrieving information from matrixindexfinal and putting it into a new dataframe
for (i in 1:(length(dades2[,1]))) {
  for (j in 1:length(z[,1])) {
    if (substr(dades2[i,1],1,7)==substr(z[j,4],1,7)) {
      dades2[i,2]<-paste(z[j,10])
      dades2[i,9]<-paste(z[j,3])
      dades2[i,10]<-paste(z[j,6])
      dades2[i,11]<-paste(z[j,7])
    }
    else {
    }
  }
}




### Removing multiple entries for multiple matrices made from submatrices of the one in the paper
# (It has a counter for when the next entry is the same and only acts if there is a difference)
j=1
p=0
k=1
simsto<-data.frame(matrix(0,nrow=length(dades2[,1]),ncol=15))
for (i in 1:length(dades2[,1])) {
  if (substr(dades2[i,1],1,7)==substr(dades2[j,1],1,7)) {
    p=p+1
    }
  else{
    simsto[k,]<-dades2[j,]
    simsto[k,1]<-substr(dades2[j,1],1,7)
    p=p+1
    j=p
    k=k+1
  }
}
simsto[k,]<-dades2[j,]
simsto[k,1]<-substr(dades2[j,1],1,7)

### Creates a new data frame which has removed all the empty rows in the simsto
redsimsto<-simsto[simsto[,1]!=0,c(1,2,3,4,6,7,8,9,10,11,12)]
colnames(redsimsto)<-c("Paper","Organism","Raw data avilable", "Gmatrix","cormatrix","directory","data link","Reference","Taxon1","Taxon2","Multiple Species")

setwd("G:/GIThub/Rproject")
write.csv(file="pardatdes.csv", x=redsimsto,row.names=FALSE)
write.csv(file="dades2.csv", x=dades2,row.names=FALSE)



############################# Checking the list of downloaded papers against the paper identifier in the pardatdesmanmod and reducing pardatdesmanmod to minimum number of rows.
pdfz<-list.files(path="G:/GIThub/Papers")
tempsto<-read.csv("C:/Users/s4284361/Documents/GitHub/Rproject/pardatdesmanmod.csv",header=TRUE,stringsAsFactors=FALSE)
grepl(substr(pdfz[1],1,7), tempsto[,1])

######Creates a new matrix where both the cov and cor column are 1
######if both matrices retrieved from one paper

tempsto2<-data.frame(matrix(0,nrow=length(tempsto[,1]),ncol=15))
for (i in 1:length(tempsto[,1])) {
  k=0
  for (j in 1:length(tempsto[,1])) {
    if (substr(tempsto[i,1],1,7)==substr(tempsto[j,1],1,7) & i!=j) {
      tempsto2[i,]<-tempsto[i,]
      tempsto2[i,c(4,5)]<-c(1,1)
      k=k+1
    }
    else {
      
    }
  }
  if (k==0) {
    tempsto2[i,]<-tempsto[i,]
  }
  else {
    
  }
}


tempsto2<-tempsto2[!duplicated(tempsto2[,1]),] ### Check which entrys have duplicates and makes a dataframe without duplicates
tempsto3<-tempsto2[order(tempsto2[,1]),] ### Orders the matrix by the paper name string
tempsto3[,16]<-rep(0,length=length(tempsto3[,1]))
for (i in 1:length(tempsto3[,1])) {
  for (j in 1:length(pdfz)) {
    if (substr(tempsto3[i,1],1,7)==substr(pdfz[j],1,7)) {
      tempsto3[i,16]<-1
    }
    else {
      
    }
  }
}
setwd("C:/Users/s4284361/Documents/GitHub/Rproject")
write.csv(file="OrdDatDes.csv", x=tempsto3,row.names=FALSE)
tempsto3[tempsto3[,16]==0,][,1] #### Gets the subset of papers without a corresponding pdf


#### Used to make sure a pdf was retrieved for every entry.
}

###### Processing Gcor and Gcov matrices from Pritcher's paper.

### First need to create list of path files. Making it a function to avoid clutter in work space
#dir1<-"C:/Users/s4284361/Documents/GitHub/Pitchers_PTRS2014/Data/Gmats_&_means_as_CSVs/"
#dir2<-"C:/Users/s4284361/Documents/GitHub/Pitchers_PTRS2014/Data/Gmats_Cor_as_CSVs/"
## Creates list of pathfiles for files in a folder.
path_file <- function(dir1,dir2) {
  if (missing(dir2)) {
    p<-paste(dir1,list.files(path=dir1),sep="/")
    return (list(dir1=p))
  }
  else {
  q<-paste(dir2,list.files(path=dir2),sep="/")
  p<-paste(dir1,list.files(path=dir1),sep="/")
  return (list(dir1=p,dir2=q))
  }
}
#### later found that full.names is an option to get complete path_file for documents in a directory for list.files/dir.
## From Pitcher's, create list of matrices for use in R modified to by functions.
MatasList<-function(dir1,dir2){
if (missing(dir2)) {
  q<-path_file(dir1)
  setwd(dir1)
  matrices <- dir()
  no.mats <- length(q[[1]][grepl("csv", q[[1]])])
  matrix_list <- list()
  
  # this loop reads in each matrix from the folder of .csv files and writes
  # them into a list
  for (i in 1:no.mats) {
    matrix_list[[i]] <- read.csv(q[[1]][grepl("csv", q[[1]])][i])
  }
  names(matrix_list) <- matrices[grepl("csv",q[[1]])]
  return(matrix_list)
}  
else {
q<-path_file(dir1,dir2)
setwd(dir1)
matrices <- dir()
no.mats <- length(q[[1]][grepl("csv", q[[1]])])
matrix_list <- list()

# this loop reads in each matrix from the folder of .csv files and writes
# them into a list
for (i in 1:no.mats) {
  matrix_list[[i]] <- read.csv(q[[1]][grepl("csv", q[[1]])][i])
}
names(matrix_list) <- matrices[grepl("csv",q[[1]])]

#### Now for correlation matrices
## From Pitcher's, create list of matrices for use in R.
setwd(dir2)
matricesc <- dir()
no.matsc <- length(q[[2]][grepl("csv", q[[2]])])
matrix_listc <- list()
# this loop reads in each matrix from the folder of .csv files and writes
# them into a list
for (i in 1:no.matsc) {
  matrix_listc[[i]] <- read.csv(q[[2]][grepl("csv", q[[2]])][i])
}
names(matrix_listc) <- matricesc[grepl("csv",q[[2]])]
return (list(matrix_listc=matrix_listc,matrix_list=matrix_list))
}
}


#### Write csv files from list function with appropiate name

WriteMatList<-function(list,dir) {
  for (i in 1:length(list)) {
    setwd(dir)
    write.csv(file=names(list[i]), x=list[[i]],row.names=FALSE)
  }
  
  
}

MaxnoTraits<- function(list) {
  q<-0
  for (i in list) {
    p<-length(i[1,])
    q<-max(p,q)
  }
  return(q)
}
#### Projection of P through G after making submatrices and standardising in some manner. ### There is an implicit assumption of equal length between the two lists
##PthroughG<-function(dir, dir2) { ## dir is the directory where writeMatList has written the new P matrices and dir2 is the dir containing the G matrices
  ##Pmat<-MatasList(dir)
  ##Gmat<-MatasList(dir2)
  ##ProjPG<-matrix(NA,nrow=length(Pmat),ncol=MaxnoTraits(Pmat))
  ##for (i in 1:length(Pmat)) {
    ##Psto<-Pmat[grepl(names(Pmat[i]),names(Gmat))]
    ##Gsto<-Gmat[grepl(names(Pmat[i]),names(Gmat))]
    ##for (j in 1:length(Psto[[1]][1,])){
      ##ProjPG[i,j]<-t((eigen(Psto[[1]])$vectors[,j]))%*%Gsto[[1]]%*%(eigen(Psto[[1]])$vectors[,j])
    ##}
  ##}
  ##return(ProjPG)
##}

#### Create list of empty csv file to later be filled with pmatrices from paper. Named files are made based on OrdDatDes pmatrix column.
#blank_csvs<-function(dir1,dir2) { #### "H:/OrdDatDes.csv"  changed directory of data due to permission issues
#  z<-read.csv(dir1,header=TRUE,stringsAsFactors=FALSE)
#  q<-z[z[,12]>0,] ### subset of lines with pmatrices
#  for (i in 1:length(q[,1])) {
#    setwd(dir2)
#    write.csv(paste(q[i,1],".1.csv",sep=""),x=NULL)
#    write.csv(paste(q[i,1],".2.csv",sep=""),x=NULL)
#  }
#}
#
#blank_csvs("H:/OrdDatDes.csv","C:/Users/s4284361/Documents/GitHub/Rproject/Pmatrices")



##### Note that con1992 might not be suitable to use due to unclear source of correlation matrix from Pitcher's data
##### Con2003 was not currently used due to Vp =Vg+Ve not matching the Vp mentioned and table, might be due to descrepancy in the degrees of freedom of each. Will attempt to address this tomorrow.
#### Del1995 omitted autocorrelations, this leaves some blanks in the matrix, it appears Pitcher's replaced these with 0s. Instead I might remove those traits as adding those zeros might affect the results.
### Kau2003 lacked information in matrix for the phenotypic correlations



#### Function for creating a matrix index for the gathered Pmatrices which will be used to create all of the relevant submatrices for each corresponding G matrix.
## This might be done in two steps, first create an index for the Pmatrices, enter the relevant G matrice manually and then extract a subset of Pitcher's matrixindex to be able to select by taxon and trait_type.

Pmatindex<-function(dir) { ## dir ="C:/Users/s4284361/Documents/GitHub/Rproject/Pmatrices"
  p<-list.files(path=dir)
  p<-p[grepl(".csv",p)] ### Takes only the csv files
  q<-read.csv("H:/OrdDatDes.csv",header=TRUE,stringsAsFactors=FALSE) ## Reads in orderdatamatrix
  z<-data.frame(matrix(NA,nrow=length(p),ncol=length(q)))
  for (i in 1:length(p)) {
    z[i,]<-q[grepl(substr(p[[i]],1,7),q[,1]),]
    z[i,1]<-substr(p[[i]],1,9)
  }
  nam<-names(q)
  nam[1]<-"Pmatrix"
  names(z)<-nam
  z$GmatId<-rep(NA,length(z[,16]))
  setwd("C:/Users/s4284361/Documents/GitHub/Rproject")
    write.csv("Pmatindex.csv",x=z,row.names=FALSE)
}

##### Function which checks GmatId, get trait_type from MatrixIndexFinal and appends it as another column to PMatindex

### Pmatindex needs to be modified a bit before this is used 8/01/15 11:18 a.m 
### Pmatindex modification made 8/01/15 11:54 a.m
# At work comp dir1<-"C:/Users/s4284361/Documents/GitHub/Rproject/Pmatindex.csv"
# at work comp dir 2<-"C:/Users/s4284361/Documents/GitHub/Pitchers_PTRS2014/Data/MatrixIndexFinal.csv"
Trait_append<-function(dir1,dir2) {
  q<-read.csv(dir1,stringsAsFactors=FALSE) ###Pmatindex
  l<-read.csv(dir2,stringsAsFactors=FALSE) ####MatrixINdexFinal 
  q$trait.type<-rep(NA,length(q[,1]))
  for (i in 1:length(q[,1])) {
    q[i,18]<-l[grepl(q[i,17],l[,5]),][8]
  }
  setwd("C:/Users/s4284361/Documents/GitHub/Rproject")
  write.csv("PmatIndex.csv",x=q,row.names=FALSE)
}
### Appends trait.type and overwrite old Pmatindex
###### New Submatrix function which check GmatId to compare Pmat to that Gmat

#### Counter function
matdimcount<-function(list) {
  p<-rep(0,length(list))
    for (i in 1:length(list)) {
      p[i]<-length(list[[i]][[1]])
    }
    q<-rep(0,max(p))
    z<-rep(0,max(p))
    for (i in 1:max(p)) {
      q[i]<-sum(p>=i)
      z[i]<-sum(p==i)
    }
    return (list(DimNum=p,DimGre=q,DimEqu=z))
}


#### Now need selecter function to extract relevant Gmatrices based on Pmatindex
MatfromInd<-function(dir1,dir2){ ## dir1 <- cor matrix folder ##dir 2 PmatIndex file_file
    q<-path_file(dir1)
    setwd(dir1)
    matrices <- dir()
    no.mats <- length(q[[1]][grepl(".csv", q[[1]])])
    matrix_list <- list()
    
    # this loop reads in each matrix from the folder of .csv files and writes
    # them into a list
    for (i in 1:no.mats) {
      matrix_list[[i]] <- read.csv(q[[1]][grepl(".csv", q[[1]])][i])
    }
    names(matrix_list) <- matrices[grepl(".csv",q[[1]])]
    z<-read.csv(dir2,stringsAsFactors=FALSE)
    GPmat<-list()
    for (i in 1:length(z[,17])) {
      GPmat[[i]]<-matrix_list[grepl(z[i,17],names(matrix_list))][[1]]
      names(GPmat)[i]<-z[i,17]
    }
    return (GPmat) ### return list of list of list
}

##### submatrix function redone to accept paired data based on PmatIndex

### Current directories dir1<-"C:/Users/s4284361/Documents/GitHub/Rproject/Pmatrices"
# dir2<-"C:/Users/s4284361/Documents/GitHub/Pitchers_PTRS2014/Data/Gmats_Cor_as_CSVs"
#dir3<-"C:/Users/s4284361/Documents/GitHub/Rproject/Pmatindex.csv"
#dir4<-"C:/Users/s4284361/Documents/GitHub/Rproject/Psubmatrices"

###### Example directories in new home machine
dir1<-"C:/GITSto/Rproject/Pmatrices"
  dir2<-"C:/GITSto/Pitchers/Pitchers_PTRS2014/Data/Gmats_Cor_as_CSVs"
  dir3<-"C:/GITSto/Rproject/Pmatindex.csv"
  dir4<-"C:/GITSto/Rproject/Psubmatrices"



######
Psubmats<-function(dir1,dir2,dir3,dir4) { ### dir1 - Pmatrices , Gmatrices, dir3-Pmatindex, dir4-storagelocation
  ### It is assumed dir1 is the larger matrix and dir 2 the list of smaller ones
  modmatsto=list()
  p<-MatasList(dir1)## Creates list of P  matrices
  q<-MatfromInd(dir2,dir3) ### creates list of relevant G matrices entreted in Pmatindex
  q<-q[!duplicated(names(q))]
  z<-read.csv(dir3,stringsAsFactors=FALSE)
  for (i in 1:length(z[,1])) {
    Pmat<-p[grepl(z[i,1],names(p))]
    Gmat<-q[grepl(z[i,17],names(q))] ## takes substring from dir1 and compare and extract matrix with matching name from dir 2
    matr<-data.frame(matrix(0,nrow=length(Pmat[[1]]),ncol=length(Gmat[[1]])))
    colnames(matr)<-names(Gmat[[1]])
    colnosto<-numeric(length(Gmat[[1]]))
    pnam<-gsub(".","",names(Pmat[[1]]),fixed=TRUE)
    gnam<-gsub(".","",names(Gmat[[1]]),fixed=TRUE)
    for (m in 1:length(Gmat[[1]])) {
      if (sum(grepl(gnam[m],pnam,ignore.case=TRUE))==1) {
        matr[,m]<-Pmat[[1]][,grepl(gnam[m],pnam,ignore.case=TRUE)]
        colnosto[m]<-match(TRUE,grepl(gnam[m],pnam,ignore.case=TRUE))
      }
      else if ((sum(grepl(gnam[m],pnam,ignore.case=TRUE))>1)) {
        colno<-match(gnam[m],pnam)
        matr[,m]<-Pmat[[1]][,colno]
        colnosto[m]<-match(gnam[m],pnam)
      }
    }
    matr<-matr[colnosto,]
    modmatsto[[names(Gmat)]]<-matr   #### name matrix with corresponding identifier from other list
    }
  for (j in 1:length(modmatsto)) {
    names(modmatsto)[j]<-paste(names(modmatsto)[j],".csv",sep="")
  }
  WriteMatList(modmatsto, dir4)
  return (modmatsto)
}
#### need a function to extract a submatrix from a larger matrix such that it makes the variables in another matrix
### for storage current grepl command grepl(colnames(mno[[2]][[1]])[1],colnames(Q[[1]]),ignore.case=TRUE)
### modified to grepl(gsub(".","",colnames(mno[[2]][[1]])[3],fixed=TRUE),colnames(Q[[1]]),ignore.case=TRUE) to counter different header style

##### Need to develop and change P through G function to accept new data
### takes Psubmatrices and turns all of them into correlation matrices
#cor<-MatasList(dir4)
#cor<-listcov2cor(cor)
#WriteMatList(cor,dir4)
###
PthroughG<-function(dir4,dir2,dir3) { ## dir1 is the directory of P submatrices and dir2 is the dir containing the G matrices
  Pmats<-MatasList(dir4)
  Gmats<-MatfromInd(dir2,dir3) ### creates list of relevant G matrices entreted in Pmatindex
  Gmats<-Gmats[!duplicated(names(Gmats))]
  Gmats<-Gmats[sort(names(Gmats))]
  ProjPG<-matrix(NA,nrow=length(Pmats),ncol=MaxnoTraits(Pmats))
  Peigsto<-matrix(NA,nrow=length(Pmats),ncol=MaxnoTraits(Pmats))
  Geigsto<-matrix(NA,nrow=length(Pmats),ncol=MaxnoTraits(Pmats))
  nullnames<-rep(0,length(Pmats))
  rownames(ProjPG)<-nullnames
  rownames(Peigsto)<-nullnames
  rownames(Geigsto)<-nullnames
  for (i in 1:length(Pmats)) {
    Psto<-Pmats[grepl(names(Gmats[i]),names(Pmats))]
    Gsto<-Gmats[grepl(names(Gmats[i]),names(Pmats))]
    rownames(ProjPG)[i]<-names(Gmats[i])
    rownames(Peigsto)[i]<-names(Gmats[i])
    rownames(Geigsto)[i]<-names(Gmats[i])
    for (j in 1:length(Psto[[1]][,1])){
      ProjPG[i,j]<-t((eigen(Psto[[1]])$vectors[,j]))%*%as.matrix(Gsto[[1]])%*%(eigen(Psto[[1]])$vectors[,j])
      Peigsto[i,j]<-eigen(Psto[[1]])$values[j]
      Geigsto[i,j]<-eigen(Gsto[[1]])$values[j]
    }
  }
  return(list(ProjPG=ProjPG,Peigsto=Peigsto,Geigsto=Geigsto))
}

### Function to manually do cov2cor since cov2cor does nto appear to work.
listcov2cor<-function(list2) {
  for (i in 1:length(list2)) {
    if (sum(diag(as.matrix(list2[[i]])))!=length(list2[[i]])) {
      varmat<-matrix(0,nrow=length(list2[[i]]),ncol=length(list2[[i]]))
      var<-diag(as.matrix(list2[[i]]))
      for (j in 1:length(list2[[i]])) {
        for (k in 1:length(list2[[i]])) {
          varmat[j,k]<-sqrt(var[j])*sqrt(var[k])
      
        }
      }
      list2[[names(list2[i])]]<-as.matrix(list2[[i]])/varmat
    }
    else {
      
    }
  }
  return (list2)
}

###### Function to extract only the data off certain length based on the no . of NAs in the row. Function is not going to be made scalable.
sizeddata<-function(list2,size) {
ProjPG<-list2[[1]]
Peigsto<-list2[[2]]
Geigsto<-list2[[3]]
SProj<-data.frame(matrix(NA,nrow=length(ProjPG[,1]),ncol=size))
SPeig<-data.frame(matrix(NA,nrow=length(ProjPG[,1]),ncol=size))
SGeig<-data.frame(matrix(NA,nrow=length(ProjPG[,1]),ncol=size))
counter<-1
  for (i in 1:length(ProjPG[,1])){
    if (is.na(match(NA,ProjPG[i,]))) {
      
    }
    else if (match(NA,ProjPG[i,])==(size+1)) {
      SProj[counter,]<-ProjPG[i,1:size]
      SPeig[counter,]<-Peigsto[i,1:size]
      SGeig[counter,]<-Geigsto[i,1:size]
      counter<-counter+1
    }
    else {
      
    }
  }
SProj<-SProj[1:(counter-1),]
SPeig<-SPeig[1:(counter-1),]
SGeig<-SGeig[1:(counter-1),]
return(list(SProj,SPeig,SGeig))
}
#### Creating an index document for Psubmatrices
##dir="C:/GITSto/Rproject/Psubmatrices"
Psubmatindex<-function(dir) { ## dir ="C:/Users/s4284361/Documents/GitHub/Rproject/Psubmatrices"
  p<-list.files(path=dir)
  p<-p[grepl(".csv",p)] ### Takes only the csv files
  q<-read.csv("C:/GITSto/Pitchers/Pitchers_PTRS2014/Data/MatrixIndexFinal.csv",header=TRUE,stringsAsFactors=FALSE) ## Reads in orderdatamatrix
  m<-read.csv("C:/GITSto/Rproject/Pmatindex.csv",header=TRUE,stringsAsFactors=FALSE)
  z<-data.frame(matrix(NA,nrow=length(p),ncol=(length(q)+2)))
  for (i in 1:length(p)) {
    z[i,]<-q[grepl(strsplit(p[[i]],".csv"),q[,5]),]
    z[i,15]<-m[grepl(strsplit(p[[i]],".csv"),m[,17]),11]
    z[i,16]<-m[grepl(strsplit(p[[i]],".csv"),m[,17]),14]
  }
  nam<-names(q)
  names(z)<-nam
  names(z)[c(15,16)]<-c("Title","DOI")
  setwd("C:/GITSto/Rproject")
  write.csv("Psubmatindex.csv",x=z,row.names=FALSE)
}

##### Extract Trait_type
##Directory Used MatasList("C:/Users/s4284361/Documents/GitHub/Rproject/Psubmatrices")

#####
### Current directories dir1<-"C:/Users/s4284361/Documents/GitHub/Rproject/Pmatrices"
# dir2<-"C:/Users/s4284361/Documents/GitHub/Pitchers_PTRS2014/Data/Gmats_Cor_as_CSVs"
#dir3<-"C:/Users/s4284361/Documents/GitHub/Rproject/Pmatindex.csv"
#dir4<-"C:/Users/s4284361/Documents/GitHub/Rproject/Psubmatrices"

#####Trialling final graph

#### P & G mats = list (will use matfromind and mataslist functions) ## Ending at uni = C:/Users/s4284361/Documents/GitHub at home =G:/GIThub
PaG<-function() {
  q<-MatasList("C:/GITSto/Rproject/Psubmatrices")
  p<-ReMatfromInd("C:/GITSto/Pitchers/Pitchers_PTRS2014/Data/Gmats_Cor_as_CSVs","C:/GITSto/Rproject/Psubmatindex.csv")
  return  (list(Pcor=q,Gcor=p))
}

####
matsubsample<-function(PG,T_no) {
  q<-PG[[1]]
  p<-PG[[2]]
  sto<-list(Psub=list(),Gsub=list())
  k=1
  for (i in 1:length(q)) {
    if (length(q[[i]][1,])<T_no) {
      
    }
    else if (length(q[[i]][1,])==T_no) {
      sto[[1]][[k]]<-q[[i]]
      names(sto[[1]])[k]<-strsplit(names(q[i]),".csv")
      sto[[2]][[k]]<-p[[i]]
      names(sto[[2]])[k]<-strsplit(names(p[i]),".csv")
      k=k+1
    }
    else if (length(q[[i]][1,])>T_no) {
      randno<-sample(1:length(q[[i]][1,]),T_no)
      sto[[1]][[k]]<-q[[i]][randno, randno]
      names(sto[[1]])[k]<-strsplit(names(q[i]),".csv")
      sto[[2]][[k]]<-p[[i]][randno, randno]
      names(sto[[2]])[k]<-strsplit(names(p[i]),".csv")
      k=k+1
    } 
    else {      
    }
  }
  return (sto)
}
#### command used for above function ReMatfromInd("C:/Users/s4284361/Documents/GitHub/Pitchers_PTRS2014/Data/Gmats_Cor_as_CSVs","C:/Users/s4284361/Documents/GitHub/Rproject/Psubmatindex.csv")

### Redone MatfromInd For Psubmatfile

ReMatfromInd<-function(dir1,dir2){ ## dir1 <- cor matrix folder ##dir 2 PmatIndex file_file
  q<-path_file(dir1)
  setwd(dir1)
  matrices <- dir()
  no.mats <- length(q[[1]][grepl(".csv", q[[1]])])
  matrix_list <- list()
  
  # this loop reads in each matrix from the folder of .csv files and writes
  # them into a list
  for (i in 1:no.mats) {
    matrix_list[[i]] <- read.csv(q[[1]][grepl(".csv", q[[1]])][i])
  }
  names(matrix_list) <- matrices[grepl(".csv",q[[1]])]
  z<-read.csv(dir2,stringsAsFactors=FALSE)
  GPmat<-list()
  for (i in 1:length(z[,5])) {
    GPmat[[i]]<-matrix_list[grepl(z[i,5],names(matrix_list))][[1]]
    names(GPmat)[i]<-z[i,5]
  }
  return (GPmat) ### return list of list of list
}

#### modified P through G to accept list instead of dir

PthroughG2<-function(tests) { ## dir1 is the directory of P submatrices and dir2 is the dir containing the G matrices
  Pmats<-tests[[1]]
  Gmats<-tests[[2]] ### creates list of relevant G matrices entreted in Pmatindex
  ProjPG<-matrix(NA,nrow=length(Pmats),ncol=MaxnoTraits(Pmats))
  Peigsto<-matrix(NA,nrow=length(Pmats),ncol=MaxnoTraits(Pmats))
  Geigsto<-matrix(NA,nrow=length(Pmats),ncol=MaxnoTraits(Pmats))
  nullnames<-rep(0,length(Pmats))
  rownames(ProjPG)<-nullnames
  rownames(Peigsto)<-nullnames
  rownames(Geigsto)<-nullnames
  for (i in 1:length(Pmats)) {
    Psto<-Pmats[grepl(names(Gmats[i]),names(Pmats))]
    Gsto<-Gmats[grepl(names(Gmats[i]),names(Pmats))]
    rownames(ProjPG)[i]<-names(Gmats[i])
    rownames(Peigsto)[i]<-names(Gmats[i])
    rownames(Geigsto)[i]<-names(Gmats[i])
    for (j in 1:length(Psto[[1]][,1])){
      ProjPG[i,j]<-t((eigen(Psto[[1]],symmetric=TRUE)$vectors[,j]))%*%as.matrix(Gsto[[1]])%*%(eigen(Psto[[1]],symmetric=TRUE)$vectors[,j])
      Peigsto[i,j]<-eigen(Psto[[1]],symmetric=TRUE)$values[j]
      Geigsto[i,j]<-eigen(Gsto[[1]],symmetric=TRUE)$values[j]
    }
  }
  return(list(ProjPG=ProjPG,Peigsto=Peigsto,Geigsto=Geigsto))
}

#### Analysis changed to ignore trait_types. Variability was too large to discriminate between groups. Also, con1992.170 was removed due to correlation values far greater than 1.
if (FALSE) {

Pcor<-listcov2cor(MatasList("G:/GIThub/Rproject/Psubmatrices"))
WriteMatList(Pcor,"G:/GIThub/Rproject/Psubmatrices")

}

MatasColumn<-function(Pmat,Gmat,Index1) { ### Takes lists of matrices and turns them into a column form
  ## Note that they must be the same length for this to work.
  ##It assumes that your are using the final form of Paired P and G matrices
  Tmax=MaxnoTraits(Pmat)
  ColumnSto<-data.frame(matrix(NA,nrow=(length(Pmat)*Tmax),ncol=9))
  colnames(ColumnSto)<-c("Matrix","P&GmatID","Organism","Population","Trait1", "Trait2","Pcorr","Gcorr","AbsDiff")
  Ind<-read.csv(Index1,header=TRUE,stringsAsFactors=FALSE)
  rowcount<-1
  for (i in 1:length(Pmat)) {
    Test<-(length(Pmat[[i]])-1)
    for (j in 1:(length(Pmat[[i]])-1)) {
      for (k in 1:Test) {
      ColumnSto[rowcount,2]<-names(Pmat[i])
      ColumnSto[rowcount,3]<-Ind[grepl(strsplit(names(Pmat[i]),".csv"),Ind[,5]),10]###Organism name from index
      ColumnSto[rowcount,5]<- names(Pmat[[i]])[j]### Trait name 1 (Column) something like names(Pmat[[i]][,j])
      ColumnSto[rowcount,6]<- names(Pmat[[i]])[k+j]## Trait name 2 (row)
      ColumnSto[rowcount,7]<- Pmat[[i]][j,k+j]
      ColumnSto[rowcount,8]<- Gmat[[i]][j,k+j]
      ColumnSto[rowcount,9]<-abs(ColumnSto[rowcount,7]-ColumnSto[rowcount,8]) ## Absolute difference
      rowcount<-rowcount+1
      }
      Test<-Test-1
    }
  }
  ColumnSto<-ColumnSto[!is.na(ColumnSto[,8]),]
  return(ColumnSto)
  
}
#### Writing ColumSto as csv file for other use.
#PG<-PaG()
ColumnSto<-MatasColumn(PG[[1]],PG[[2]],"C:/GITSto/Rproject/Psubmatindex.csv")
write.csv("MatricesAsColumn.csv",x=ColumnSto,row.names=FALSE)
###
PG<-PaG() ##### note gem2001.473 was earlier removed. Rerunning everything includes it. Must be removed for later functions to work. 
tests<-matsubsample(PG,5)
test2<-PthroughG2(tests)
if (FALSE) {
lkt<-tests
names(lkt[[1]])<-paste(names(lkt[[1]]),".csv",sep="")
names(lkt[[2]])<-paste(names(lkt[[2]]),".csv",sep="")
WriteMatList(lkt[[1]],"C:/Users/s4284361/Documents/GitHub/Rproject/Psubsampledmatrices")
WriteMatList(lkt[[2]],"C:/Users/s4284361/Documents/GitHub/Rproject/Gsubsampledmatrices")


Gdir<-"C:/Users/s4284361/Documents/GitHub/Rproject/Graphs"

mypath<-file.path(paste(Gdir),paste("All3","_boxplot" ,paste("T5n=43"), ".pdf", sep = ""))
pdf(file=mypath)
boxplot(test2[[2]],col="grey",boxwex=0.18,xaxis=NULL,xlab="Eigenvector",ylab="Scaled variance",ylim=c(min(test2[[2]],test2[[3]],test2[[1]]),max(test2[[2]],test2[[3]],test2[[1]])))
boxplot(test2[[3]],add=TRUE, at=1.2:(5+0.2), boxwex=0.18,xaxt='n')
boxplot(test2[[1]],col="gray42",add=TRUE, at=1.4:(5+0.4), boxwex=0.18,xaxt='n')
legend(x="topright", c("Phenotypic_Correlation","Genetic_Correlation","Projection of P through G"), fill=c("grey","white","gray42"))
dev.off()

z<-(test2$Peigsto-test2$Geigsto)

mypath<-file.path(paste(Gdir),paste("difP&G","_boxplot" ,paste("T5n=43"), ".pdf", sep = ""))
pdf(file=mypath)
boxplot(z,col="grey",boxwex=0.18,xaxis=NULL,xlab="Trait or Eigenvector",ylab="Difference in scaled variance")
legend(x="topleft", c("eig(P)-eig(G)"), fill=c("grey"))
abline(0,0,lty=2)
dev.off()


##Differences between P and G for Trait types
mn<-data.frame(matrix(0,nrow=43,ncol=16))
z<-read.csv("C:/Users/s4284361/Documents/GitHub/Rproject/Psubmatindex.csv",header=TRUE,stringsAsFactors=FALSE)
for (i in 1:length(names(tests[[1]]))) {mn[i,]<-z[grepl(names(tests[[1]])[i],z$filename),]}

PsubtraitM<-tests[[1]][grepl("M",mn[,8])]
PsubtraitS<-tests[[1]][grepl("S",mn[,8])]
PsubtraitL<-tests[[1]][grepl("L",mn[,8])]
GsubtraitM<-tests[[2]][grepl("M",mn[,8])]
GsubtraitS<-tests[[2]][grepl("S",mn[,8])]
GsubtraitL<-tests[[2]][grepl("L",mn[,8])]

subM<-list(PsubtraitM,GsubtraitM)
subS<-list(PsubtraitS,GsubtraitS)
subL<-list(PsubtraitL,GsubtraitL)

MPtG<-PthroughG2(subM)
SPtG<-PthroughG2(subS)
LPtG<-PthroughG2(subL)

DifM<-(MPtG[[2]]-MPtG[[3]])
DifS<-(SPtG[[2]]-SPtG[[3]])
DifL<-(LPtG[[2]]-LPtG[[3]])

mypath<-file.path(paste(Gdir),paste("SplitT_P","_boxplot" ,paste("T5n=43"), ".pdf", sep = ""))
pdf(file=mypath)
boxplot(DifL,col="white",boxwex=0.18,xaxis=NULL,xlab="Eigenvector",ylab="eig(P) - eig(G)",ylim=c(min(DifM,DifL,DifS),max(DifM,DifL,DifS)))
boxplot(DifM,col="grey",add=TRUE, at=1.2:(5+0.2), boxwex=0.18,xaxt='n')
boxplot(DifS,col="gray42",add=TRUE, at=1.4:(5+0.4), boxwex=0.18,xaxt='n')
abline(0,0,lty=2)
text(x=c(0.5,0.5),y=c(0.25,-0.25), labels=c("P>G","P<G"))
legend(x="bottomright", c("Life History - 7","Morphology - 26","Sexually Selected - 10"), fill=c("white","grey","gray42"))
dev.off()

### INcluding compiled dataset
z<-(test2$Peigsto-test2$Geigsto)

mypath<-file.path(paste(Gdir),paste("SplitT_P&All","_boxplot" ,paste("T5n=43"), ".pdf", sep = ""))
pdf(file=mypath)
boxplot(DifL,col="white",boxwex=0.18,xaxis=NULL,xlab="Trait or Eigenvector",ylab="eig(P) - eig(G)",ylim=c(-1.7,1))
boxplot(DifM,col="grey",add=TRUE, at=1.2:(5+0.2), boxwex=0.18,xaxt='n')
boxplot(DifS,col="gray42",add=TRUE, at=1.4:(5+0.4), boxwex=0.18,xaxt='n')
boxplot(z,col="black",add=TRUE, at=1.6:(5+0.6), boxwex=0.18, xaxt='n')
abline(0,0,lty=2)
legend(x="bottomright", c("Life History - 7","Morphology - 26","Sexually Selected - 10", "All -43"), fill=c("white","grey","gray42","black"))
dev.off()

mypath<-file.path(paste(Gdir),paste("SplitT_PtG","_boxplot" ,paste("T5n=43"), ".pdf", sep = ""))
pdf(file=mypath)
boxplot(LPtG[[1]],col="white",boxwex=0.18,xaxis=NULL,xlab="Trait or Eigenvector",ylab="P'GP",ylim=c(min(MPtG[[1]][,5]),max(MPtG[[1]][,1])))
boxplot(MPtG[[1]],col="grey", add=TRUE, at=1.2:(5+0.2), boxwex=0.18,xaxt='n')
boxplot(SPtG[[1]],col="gray42",add=TRUE, at=1.4:(5+0.4), boxwex=0.18,xaxt='n')
legend(x="topright", c("Life History - 7","Morphology - 26","Sexually Selected - 10"), fill=c("white","grey","gray42"))
abline(1,0,lty=2)
dev.off()


#### Splitting data by plant/animal taxon
PtaxonA<-tests[[1]][grepl("A",mn[,7])]
PtaxonP<-tests[[1]][grepl("P",mn[,7])]
GtaxonA<-tests[[2]][grepl("A",mn[,7])]
GtaxonP<-tests[[2]][grepl("P",mn[,7])]

subA<-list(PtaxonA,GtaxonA)
subP<-list(PtaxonP,GtaxonP)

APtG<-PthroughG2(subA)
PPtG<-PthroughG2(subP)

DifA<-(APtG[[2]]-APtG[[3]])
DifP<-(PPtG[[2]]-PPtG[[3]])

mypath<-file.path(paste(Gdir),paste("SplitTaxon","_boxplot" ,paste("T5n=43"), ".pdf", sep = ""))
pdf(file=mypath)
boxplot(DifA,col="white",boxwex=0.18,xaxis=NULL,xlab="Trait or Eigenvector",ylab="eig(P) - eig(G)")
boxplot(DifP,col="grey",add=TRUE, at=1.2:(5+0.2), boxwex=0.18,xaxt='n')
abline(0,0,lty=2)
legend(x="bottomright", c("Animal -32","Plant -11"), fill=c("white","grey"))
dev.off()

mypath<-file.path(paste(Gdir),paste("SplitTaxonPGP","_boxplot" ,paste("T5n=43"), ".pdf", sep = ""))
pdf(file=mypath)
boxplot(APtG[[1]],col="white",boxwex=0.18,xaxis=NULL,xlab="Trait or Eigenvector",ylab="P'GP")
boxplot(PPtG[[2]],col="grey",add=TRUE, at=1.2:(5+0.2), boxwex=0.18,xaxt='n')
abline(0,0,lty=2)
legend(x="topright", c("Animal -32","Plant -11"), fill=c("white","grey"))
dev.off()


nullmat<-matrix(0,nrow=length(test2[[1]][,1]),ncol=15)
nullmat[,1:5]<-test2[[1]]
nullmat[,6:10]<-test2[[2]]
nullmat[,11:15]<-test2[[3]]

#### Unused function and graph
anglePG<-function(lists){
  angle<-matrix(0,nrow=length(lists[[1]]),ncol=5*length(lists[[1]][[1]][1,]))
  for (i in 1:length(lists[[1]])) {
    for (j in 1:length(lists[[1]][[1]][1,])){
      for (k in 1:5) {
    angle[i,(j-1)*5+k]<-acos(t(eigen(lists[[1]][[i]])$vectors[,j])%*%eigen(lists[[2]][[i]])$vectors[,k]/(sqrt(sum(eigen(lists[[1]][[i]])$vectors[,j]*eigen(lists[[1]][[i]])$vectors[,j])*sum(eigen(lists[[2]][[i]])$vectors[,k]*eigen(lists[[2]][[i]])$vectors[,k]))))*180/pi  
      }
    }
  }
  return(angle)
}
q<-angle

boxplot(q[,1:5],boxwex=0.12,col="grey",at=0.8:4.8,ylab="angle between eigenvectors of P and G",xlab="eigenvector of P")
boxplot(q[,6:10],boxwex=0.12,at=1:5-0.05,add=TRUE,col="white",xaxt='n')
boxplot(q[,11:15],add=TRUE,boxwex=0.12,at=1:5+0.10,col="grey42",xaxt='n')
boxplot(q[,16:20],add=TRUE,boxwex=0.12,at=1:5+0.25,col="tan",xaxt='n')
boxplot(q[,21:25],add=TRUE,boxwex=0.12,at=1:5+0.40,col="azure",xaxt='n')
abline(90,0,lty=2) ### angle for orthogonal matrices
cov2cor(var(q[,c(1,6,11,16,21)])) # angle correlation


}
#####


## Quick way to test if all names existi n matindex final 
#for (i in 1:length(p)) {
 # z[i,]<-sum(grepl(strsplit(p[[i]],".csv"),q[,5]))
  #z[i,15]<-sum(grepl(strsplit(p[[i]],".csv"),m[,17]))
  #z[i,16]<-sum(grepl(strsplit(p[[i]],".csv"),m[,17]))
#}

#sum(z[,1])/length(z[,1]) #!=1 indicates mistake
