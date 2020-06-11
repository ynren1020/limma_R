##################2019-08-27###########################
##check novel lncRNAs with SE is DEGs or not ##########
#######################################################

library(tidyverse)

input3<-"VCaP_SE_novel4specific_20190827.txt"
input2<-"prad.meta.vs.normal.degs.normalfiltered.loose.txt"

enhancer<-read.delim(input3,header=FALSE,stringsAsFactors = FALSE)
degs<-read.delim(input2,header=TRUE,stringsAsFactors = FALSE)

degs.sub<-degs[grep("G.",row.names(degs),ignore.case=FALSE,fixed=TRUE),]

degs.sub.join<-degs.sub[degs.sub$geneid%in%enhancer$V1,]

##logFC > 0 adj.pval<=0.05##
degs.sub.join.sub<-degs.sub.join[degs.sub.join$logFC>1&degs.sub.join$adj.P.Val<=0.05,]

write.table(degs.sub.join.sub,"prad.meta.vs.normal.degs.normalfiltered.loose.enhancer_novel.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")

write.table(degs.sub.join.sub$geneid,"prad.meta.vs.normal.degs.normalfiltered.loose.enhancer.novelname.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep = "\t")


##super enhancer, DEGs, FPKM (genes) of TCGA##
tcga<-read.delim("tcga.enhancer.DEGs.FPKM.txt",header=FALSE,stringsAsFactors = FALSE,check.names = FALSE)
cohort<-tcga[1,1:33]
tcga.sub<-tcga[-1,-seq(3,ncol(tcga),2)]
colnames(tcga.sub)<-c("genename",unlist(cohort[1,]))
rownames(tcga.sub)<-tcga.sub$genename
tcga.sub$genename<-NULL
##Character to numeric##
tcga.sub[] <- lapply(tcga.sub, function(x) as.numeric(as.character(x)))

tcga.subT<-as.data.frame(t(tcga.sub))

##specificity##
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
#Minimum 2 tissues
fTau <- function(x)
{
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(max(x)!=0)
      {
        x <- (1-(x/max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res/(length(x)-1)
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    } 
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  } 
  return(res)
}
###***###***###

###+++###
#Z-score
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fZ <- function(x)
{	
  res<-(x-mean(x))/sd(x)
  return(res)
}
###***###***###	


tau<-as.data.frame(apply(tcga.subT,2,fTau))
zscore<-as.data.frame(apply(tcga.subT,2,fZ))
zscore$cohort<-rownames(zscore)
zscore.prad<-filter(zscore,cohort=="PRAD")

write.table(tau,"TCGA.superenhancer.DEGs.tau.txt",quote = FALSE,row.names = TRUE,col.names = FALSE,sep="\t")
write.table(zscore,"TCGA.superenhancer.DEGs.zscore.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep="\t")
write.table(zscore.prad,"TCGA.superenhancer.DEGs.zscore.prad.txt",quote = FALSE,row.names = FALSE,col.names = TRUE,sep="\t")


##srpm G.44290, G.60406, G.72529 ##
g44290<-read.delim("tcga_G.44290_junction.txt",header = TRUE,stringsAsFactors = FALSE)
g60406<-read.delim("tcga_G.60406_junction.txt",header = TRUE,stringsAsFactors = FALSE)
g72529<-read.delim("tcga_G.72529_junction.txt",header = TRUE,stringsAsFactors = FALSE)
input<-"TCGA.metainfo.txt"
meta<-read.delim(input,header = TRUE,stringsAsFactors = FALSE,sep = "\t")
##join with meta##
g44290.join<-left_join(g44290,meta,by=c("sample"="FILE_ID"))
g44290.join<-na.omit(g44290.join)
g44290.join.sub<-g44290.join[,c("PROJECT","junctions","srpm")]
g44290.join.sub$srpm<-as.numeric(g44290.join.sub$srpm)
g44290.join.sub<-g44290.join.sub%>%group_by(PROJECT)%>%summarise(srpmSum=sum(srpm))
g44290.join.sub$count<-c(1109,169,500,552)
g44290.join.sub$srpmMean<-g44290.join.sub$srpmSum/g44290.join.sub$count

g44290.join.sub2<-data.frame(PROJECT=setdiff(meta$PROJECT,g44290.join.sub$PROJECT),srpmSum=0,count=100,srpmMean=0)

g44290.join.sub.join<-bind_rows(g44290.join.sub,g44290.join.sub2)

g44290.join.sub.join$tau<-fTau(g44290.join.sub.join$srpmMean)
g44290.join.sub.join$zscore<-fZ(g44290.join.sub.join$srpmMean)

##mean obtained by sum/allsamplesize##
##nofilter, but use all sample size to average##
meta.tumor<-filter(meta,TUMOR_TYPE!="Solid Tissue Normal")
meta.tumor.sum<-meta.tumor%>%group_by(PROJECT)%>%
  summarise(count=n())

g60406.join<-left_join(g60406,meta,by=c("sample"="FILE_ID"))
g60406.join<-na.omit(g60406.join)
g60406.join.sub<-g60406.join[,c("PROJECT","srpm")]
g60406.join.sub$srpm<-as.numeric(g60406.join.sub$srpm)
g60406.join.subsum<-g60406.join.sub%>%group_by(PROJECT)%>%summarise(srpmSum=sum(srpm))

g60406.join<-full_join(g60406.join.subsum,meta.tumor.sum,by="PROJECT")
g60406.join[is.na(g60406.join[])]<-0
g60406.join$srpmMean<-g60406.join$srpmSum/g60406.join$count
##specificity##
g60406.join$tau<-fTau(g60406.join$srpmMean)
g60406.join$zscore<-fZ(g60406.join$srpmMean)

##72529##
g72529.join<-left_join(g72529,meta,by=c("sample"="FILE_ID"))
g72529.join<-na.omit(g72529.join)
g72529.join.sub<-g72529.join[,c("PROJECT","srpm")]
g72529.join.sub$srpm<-as.numeric(g72529.join.sub$srpm)
g72529.join.sub<-g72529.join.sub%>%group_by(PROJECT)%>%summarise(srpmSum=sum(srpm))


g72529.join<-full_join(g72529.join.sub,meta.tumor.sum,by="PROJECT")
g72529.join[is.na(g72529.join[])]<-0
g72529.join$srpmMean<-g72529.join$srpmSum/g72529.join$count
##specificity##
g72529.join$tau<-fTau(g72529.join$srpmMean)
g72529.join$zscore<-fZ(g72529.join$srpmMean) #5.07 PRAD however only one patient




