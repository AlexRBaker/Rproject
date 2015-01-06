
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


###### Processing Gcor and Gcov matrices from Pritcher's paper.

### First need to create list of path files. Making it a function to avoid clutter in work space
dir1<-"C:/Users/s4284361/Documents/GitHub/Pitchers_PTRS2014/Data/Gmats_&_means_as_CSVs/"
dir2<-"C:/Users/s4284361/Documents/GitHub/Pitchers_PTRS2014/Data/Gmats_Cor_as_CSVs/"
## Creates list of pathfiles for files ina folder.
path_file <- function(dir1,dir2) {
  if (missing(dir2)) {
    p<-paste(dir1,list.files(path=dir1),sep="")
    return (list(dir1=p))
  }
  else {
  q<-paste(dir2,list.files(path=dir2),sep="")
  p<-paste(dir1,list.files(path=dir1),sep="")
  return (list(dir1=p,dir2=q))
  }
}
#### later found that full.names is an option to get complete path_file for documents in a directory for list.files/dir.
## From Pitcher's, create list of matrices for use in R modified to be  function.
MatasList<-function(dir1,dir2){
if (missing(dir2)) {
  q<-path_file(dir1)
  setwd(dir1)
  matrices <- dir()
  no.mats <- length(q[[1]][!grepl("list", q[[1]])])
  matrix_list <- list()
  
  # this loop reads in each matrix from the folder of .csv files and writes
  # them into a list
  for (i in 1:no.mats) {
    matrix_list[[i]] <- read.csv(q[[1]][!grepl("list", q[[1]])][i])
  }
  names(matrix_list) <- matrices[!grepl("list",q[[1]])]
  return(matrix_list)
}  
else {
q<-path_file(dir1,dir2)
setwd(dir1)
matrices <- dir()
no.mats <- length(q[[1]][!grepl("list", q[[1]])])
matrix_list <- list()

# this loop reads in each matrix from the folder of .csv files and writes
# them into a list
for (i in 1:no.mats) {
  matrix_list[[i]] <- read.csv(q[[1]][!grepl("list", q[[1]])][i])
}
names(matrix_list) <- matrices[!grepl("list",q[[1]])]

#### Now for correlation matrices
## From Pitcher's, create list of matrices for use in R.
setwd(dir2)
matricesc <- dir()
no.matsc <- length(q[[2]][!grepl("list", q[[2]])])
matrix_listc <- list()
# this loop reads in each matrix from the folder of .csv files and writes
# them into a list
for (i in 1:no.matsc) {
  matrix_listc[[i]] <- read.csv(q[[2]][!grepl("list", q[[2]])][i])
}
names(matrix_listc) <- matricesc[!grepl("list",q[[2]])]
return (list(matrix_listc=matrix_listc,matrix_list=matrix_list))
}
}

#### need a function to extract a submatrix from a larger matrix such that it makes the variables in another matrix
Psubmats<-function(dir1,dir2) {
  ### It is assumed dir1 is the larger matrix and dir 2 the list of smaller ones
  modmatsto=list()
  count=0
  p<-MatasList(dir1,dir2)
  for (c in p[[1]]) {
    Q<-p[[2]][grepl(substr(names(c),1,7),names(p[[2]]))] ## takes substring from dir1 and compare and extract matrice with matching name frmo dir 2
    for (z in Q) {
      matr<-##blank matrix
      for (i in z) {
        matr[,i]<-c[,grepl(colnames(z)[i],gsub(".","",colnames(c),fixed=TRUE),ignore.case=TRUE)]
        colnames(matr)[i]<-colnames(Q[[2]])[1] ## name traits in matrix
      }
      modmatsto[[paste(names(z))]]<-matr   #### name matrix with corresponding identifier from other list
    }
  }
  return (modmatsto)
}
### for storage current grepl command grepl(colnames(mno[[2]][[1]])[1],colnames(Q[[1]]),ignore.case=TRUE)
### modified to grepl(gsub(".","",colnames(mno[[2]][[1]])[3],fixed=TRUE),colnames(Q[[1]]),ignore.case=TRUE) to counter different header style
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
PthroughG<-function(dir, dir2) { ## dir is the directory where writeMatList has written the new P matrices and dir2 is the dir containing the G matrices
  Pmat<-MatasList(dir)
  Gmat<-MatasList(Gmat)
  ProjPG<-matrix(NA,nrow=length(Pmat),ncol=MaxnoTraits(Pmat))
  for (i in 1:length(Pmat)) {
    Psto<-Pmat[grepl(names(Pmat[i]),names(Gmat))]
    Gsto<-Gmat[grepl(names(Pmat[i]),names(Gmat))]
    for (j in 1:length(Psto[[1]][1,])){
      ProjPG[i,j]<-t((eigen(Psto[[1]])$vectors[,j]))%*%Gsto[[1]]%*%(eigen(Psto[[1]])$vectors[,j])
    }
  }
  return(ProjPG)
}

#### Might be able to include taxons if I refer to the Matrix Index, and then create two subsets of the G mat and Pmat which match the selection.
#### Taxon1 and 2 are as names in Pitcher's et al. and list1 and 2 are the list of P and G matrices to be selected from. The results should be a list of two list of matrices of equal length.
ExtractTaxon<-function(taxon1,taxon2,list1,list2,trait_type) {
  z<-read.csv("C:/Users/s4284361/Documents/GitHub/Pitchers_PTRS2014/Data/MatrixIndexFinal.csv",header=TRUE)
  #### Various decision trees for the grepl commands to get an appropiate subset of MatrixIndex to work on.
  if (missing(trait_type)) {
    
  }
  else if (missing(taxon2)) {
    
  }
  else if (missing(taxon1)) {
    
  }
  else if (missing(taxon1) & missing(taxon2)) {
    
  }
  else if (missing(taxon1) & missing(trait_type)) {
    
  }
  else if (missing(taxon2) & missing(trait_type) {
    
  }
  else if (!(missing(taxon2)|missing(taxon1)|missing(trait_type))) {
    
  }
  else {
    print("You have not entered any selection criteria")
    break
  }
  if (){
    
  }
  else {
    
  }
  ## Some output of if else tree used to get subset of p and G
  ##### Its got to take taxon 1 and 2, get the subset of the data from Z and then compare names against those in list 1 and 2 to extract the relevant matrices
  return(modlist1,modlist2)
}