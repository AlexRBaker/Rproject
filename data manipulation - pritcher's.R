p<-list.files(path="C:/Users/s4284361/Documents/GitHub/Pitchers_PTRS2014/Data/Gmats_&_means_as_CSVs")
q<-list.files(path="C:/Users/s4284361/Documents/GitHub/Pitchers_PTRS2014/Data/Gmats_Cor_as_CSVs")

datdes<-data.frame(matrix(NA,nrow=(length(p)+length(q)),ncol=9))
colnames(datdes)<-c("Paper","Organism","Raw data avilable", "Gmatrix","Pmatrix","cormatrix","directory","data link","means")
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

write.csv(file="pardatdes", x=dades2,row.names=FALSE)

setwd("C:/ABakeSumProj/Lots of data")

j=1
p=0
k=1
simsto<-data.frame(matrix(0,nrow=length(dades2[,1]),ncol=9))
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

redsimsto<-simsto[simsto[,1]!=0,c(1,2,3,4,6,7,8)]
colnames(redsimsto)<-c("Paper","Organism","Raw data avilable", "Gmatrix","cormatrix","directory","data link")
write.csv(file="pardatdes", x=redsimsto,row.names=FALSE)

setwd("C:/ABakeSumProj/Lots of data")

