
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
tempsto<-read.csv("G:/GIThub/Rproject/pardatdesmanmod.csv",header=TRUE,stringsAsFactors=FALSE)
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

tempsto3[tempsto3[,16]==0,][,1] #### Gets the subset of papers without a corresponding pdf


#### Used to make sure a pdf was retrieved for every entry.
