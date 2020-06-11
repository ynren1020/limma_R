##########2019-02-06#############
##choose genes that DE###########
##for rgt-TDF analysis###########


##k562##
dat<-read.delim("parental.degs.txt")
dat$geneid<-row.names(dat)
dat<-dplyr::filter(dat,abs(logFC)>=2&adj.P.Val<0.05)
for (i in 1:nrow(dat)){
k562.genelist[i]<-substring(dat$geneid[i],1,15)
}
write.table(k562.genelist,"k562.genelist.txt",sep="\t",col.names = FALSE,row.names = FALSE,quote = FALSE)

##vcap##
##knock down##si025
dat<-read.delim("si025.degs.txt")
dat$geneid<-row.names(dat)
dat<-dplyr::filter(dat,adj.P.Val<0.05)
#dat<-dplyr::filter(dat,abs(logFC)>=1)
#vcap<-dat$geneid
##knock down si033##
dat2<-read.delim("si033.degs.txt")
dat2$geneid<-row.names(dat2)
dat2<-dplyr::filter(dat2,adj.P.Val<0.05)
#dat<-dplyr::filter(dat,abs(logFC)>=1)
#vcap2<-dat2$geneid
dat3<-dplyr::full_join(dat,dat2,by="geneid")
dat4<-na.omit(dat3)

vcap3<-dat3$geneid
vcap4<-dat4$geneid


write.table(vcap3,"vcap3.genelist.txt",sep="\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
write.table(vcap4,"vcap4.genelist.txt",sep="\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
