p<-list.files(path="G:/GIThub/Pitchers_PTRS2014/Data/Gmats_&_means_as_CSVs")
q<-list.files(path="G:/GIThub/Pitchers_PTRS2014/Data/Gmats_Cor_as_CSVs")
z<-read.csv("G:/GIThub/Pitchers_PTRS2014/Data/MatrixIndexFinal.csv",header=TRUE)
datdes<-data.frame(matrix(NA,nrow=(length(p)+length(q)),ncol=15))
colnames(datdes)<-c("Paper","Organism","Raw data avilable", "Gmatrix","Pmatrix","cormatrix","directory","data link","Reference","Taxon1","Taxon2")


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

substr(dades2[1,1],1,7)


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

redsimsto<-simsto[simsto[,1]!=0,c(1,2,3,4,6,7,8,9,10,11,12)]
colnames(redsimsto)<-c("Paper","Organism","Raw data avilable", "Gmatrix","cormatrix","directory","data link","Reference","Taxon1","Taxon2","Multiple Species")

setwd("G:/GIThub/Rproject")
write.csv(file="pardatdes.csv", x=redsimsto,row.names=FALSE)
write.csv(file="dades2.csv", x=dades2,row.names=FALSE)
