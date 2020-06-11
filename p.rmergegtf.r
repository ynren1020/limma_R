##########Nov 28 2018##############
#######try subset cp_parental_new.gtf and cp_release_new.gtf with only gene names .p and .r####
#######then use gffcompare to get unique names for those genes#########

##load packages##
library(dplyr)
library(tidyr)

pnew<-read.delim("cp_parental_new.gtf",header = FALSE)
rnew<-read.delim("cp_release_new.gtf",header = FALSE)

pnew<-separate(pnew,V9,c("gene_id","t_id","gene_name"),sep = ";")
rnew<-separate(rnew,V9,c("gene_id","t_id","gene_name"),sep = ";")

pnewsub<-pnew[grep(".p",pnew$gene_id),]
rnewsub<-rnew[grep(".r",rnew$gene_id),]

for (i in 1:nrow(pnewsub)){
pnewsub$t_id[i]<-paste0("transcript_id ",paste(strsplit(strsplit(pnewsub$t_id[i]," ")[[1]][2],".",fixed=TRUE)[[1]][1],strsplit(strsplit(pnewsub$t_id[i]," ")[[1]][2],".",fixed=TRUE)[[1]][2],strsplit(strsplit(pnewsub$t_id[i]," ")[[1]][2],".",fixed=TRUE)[[1]][4],strsplit(strsplit(pnewsub$t_id[i]," ")[[1]][2],".",fixed=TRUE)[[1]][3],sep="."))

}

for (i in 1:nrow(pnewsub)) {

    pnewsub$gene_id[i]<-paste(strsplit(pnewsub$gene_id[i]," ")[[1]][2],paste(dQuote(strsplit(pnewsub$gene_id[i]," ")[[1]][3]),sep=""),sep=" ")
    pnewsub$t_id[i]<-paste(strsplit(pnewsub$t_id[i]," ")[[1]][1],paste(dQuote(strsplit(pnewsub$t_id[i]," ")[[1]][2]),sep=""),sep=" ")
    pnewsub$gene_name[i]<-paste(strsplit(pnewsub$gene_name[i]," ")[[1]][2],paste(dQuote(strsplit(pnewsub$gene_name[i]," ")[[1]][3]),sep=""),sep=" ")
    
}

#####need work on it###Nov28.2018
for (i in 1:nrow(pnewsub)){
    if (strsplit(pnewsub$gene_name[i]," ")[[1]][1]!="exon_number")
pnewsub$gene_name[i]<-""

}

pnewsub<-unite(pnewsub,V9,c("gene_id","t_id","gene_name"),sep = ";")
write.table(pnewsub,file="cp_parental_new_sub.gtf",row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")

#################################
for (i in 1:nrow(rnewsub)){
    rnewsub$t_id[i]<-paste0("transcript_id ",paste(strsplit(strsplit(rnewsub$t_id[i]," ")[[1]][2],".",fixed=TRUE)[[1]][1],strsplit(strsplit(rnewsub$t_id[i]," ")[[1]][2],".",fixed=TRUE)[[1]][2],strsplit(strsplit(rnewsub$t_id[i]," ")[[1]][2],".",fixed=TRUE)[[1]][4],strsplit(strsplit(rnewsub$t_id[i]," ")[[1]][2],".",fixed=TRUE)[[1]][3],sep="."))
    
}


for (i in 1:nrow(rnewsub)) {
    
    rnewsub$gene_id[i]<-paste(strsplit(rnewsub$gene_id[i]," ")[[1]][2],paste(dQuote(strsplit(rnewsub$gene_id[i]," ")[[1]][3]),sep=""),sep=" ")
    rnewsub$t_id[i]<-paste(strsplit(rnewsub$t_id[i]," ")[[1]][1],paste(dQuote(strsplit(rnewsub$t_id[i]," ")[[1]][2]),sep=""),sep=" ")
    rnewsub$gene_name[i]<-paste(strsplit(rnewsub$gene_name[i]," ")[[1]][2],paste(dQuote(strsplit(rnewsub$gene_name[i]," ")[[1]][3]),sep=""),sep=" ")
    
}

for (i in 1:nrow(rnewsub)){
    if (strsplit(rnewsub$gene_name[i]," ")[[1]][1]!="exon_number")
        rnewsub$gene_name[i]<-""
    
}
rnewsub<-unite(rnewsub,V9,c("gene_id","t_id","gene_name"),sep = ";")
write.table(rnewsub,file="cp_release_new_sub.gtf",row.names = FALSE,col.names = FALSE,quote=FALSE,sep="\t")


