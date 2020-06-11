#############################Date:2018-12-31#####################################################
#############################RNA-seq analysis-DEG#########################################
#######Input file: *_id_counts.txt files from featureCounts commands in Linux####################
#################################################################################################

####install packages "edgeR" and "limma" and "R.utils"###############
#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
#biocLite("edgeR")
####read data in and merge data together############################


####load libraries####
library("edgeR")
#library("sva")
#library("R.utils")
library("limma")
library(dplyr)
library(tidyr)
library(stringr)


# args <- commandArgs(trailingOnly = TRUE);
# input.matrix = args[1];

#nilo <- read.delim("nilo_id_counts.txt", row.names="Geneid", check.names=FALSE)
#nilosub<-nilo[,c(5,6,7)]
#nilosub<-rename(nilosub,nilo1="/home/ryang/project/k562_drug/nilo1_IN_sorted.bam")
#nilosub<-rename(nilosub,nilo2="/home/ryang/project/k562_drug/nilo2_IN_sorted.bam")

##FUNCTION TO READ SAMPLE DATA IN##
sample<-function(file){
    temp<-read.delim(file,row.names="Geneid",check.names=FALSE,skip = 1)
    tempsub<-temp[,c(5,6,7)]
    colnames(tempsub)[2]<-strsplit(strsplit(colnames(tempsub)[2],"/")[[1]][6],"_")[[1]][1]
    colnames(tempsub)[3]<-strsplit(strsplit(colnames(tempsub)[3],"/")[[1]][6],"_")[[1]][1]
    return(tempsub)
    }
nilo<-sample("nilo_id_counts_test.txt")
parental<-sample("parental_id_counts_test.txt")
release<-sample("release_id_counts_test.txt")

all<-merge(nilo,parental,by="row.names")
row.names(all)<-all$Row.names
all$Row.names<-NULL
all3<-merge(all,release,by="row.names")
row.names(all3)<-all3$Row.names
all3$Row.names<-NULL
##final merged data##
all3$Length.x<-NULL
all3$Length.y<-NULL
all3<-all3[,c(5,1,2,3,4,6,7)]
##data for DEG analysis##
data<- all3[,seq(2,7)];
group <- factor(c(1,1,2,2,3,3))


y <- DGEList(counts=data, group=group);
y$genes <- data.frame(Length=all3$Length);
rownames(y$genes) <- rownames(y$counts);

#y$genes$Symbol <- data.frame(Length=src.data$Symbol);
#rownames(y$symbol) <- rownames(y$counts);

# filtering out low expressed genes
#  the minimum number of samples in each group is two, over here.
keep <- rowSums(cpm(y)>1) >= 2
table(keep);
#keep
#FALSE  TRUE 
#44010 16765 
y <- y[keep, keep.lib.sizes=FALSE];



design <- model.matrix(~0+group);
rownames(design) <- colnames(y)
design
#group1 group2 group3
#nilo1          1      0      0
#nilo2          1      0      0
#parental1      0      1      0
#parental2      0      1      0
#release1       0      0      1
#release2       0      0      1
#attr(,"assign")
#[1] 1 1 1
#attr(,"contrasts")
#attr(,"contrasts")$group
#[1] "contr.treatment"
# Note that the filtering does not use knowledge of what treatment corresponds to each sample, so
# the filtering does not bias the subsequent differential expression analysis.
# The TMM normalization is applied to account for the compositional biases:

# normalization by the library sizes
y <- calcNormFactors(y);
y$samples;

write.table(rpkm(y), file=paste0('K562rpkm','.count.txt'), sep='\t', row.names = TRUE, quote = FALSE);

barplot(y$samples$lib.size,names=colnames(y),las=2)
# Add a title to the plot
title("Barplot of library sizes")

# normalise the read counts with 'voom' function
v <- voom(y,design,plot = TRUE)
# extract the normalised read counts
counts.voom <- v$E

boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")

# save normalised expression data into output dir
write.table(counts.voom,file="K562_counts.voom.txt",row.names=T,quote=F,sep="\t");
########################################################################################################################
# fit linear model for each gene given a series of libraries
fit <- lmFit(v, design)
##parental vs nilo##
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
matrix.2vs1 <- makeContrasts(group2-group1,levels=design)
# compute estimated coefficients and standard errors for a given set of contrasts
fit.2vs1 <- contrasts.fit(fit, matrix.2vs1)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit.2vs1 <- eBayes(fit.2vs1)
options(digits=3)

#colnames(fit$coefficients)

summary(decideTests(fit.2vs1, p.value=0.05,lfc=1))
#group2 - group1
#Down              3255
#NotSig           12472
#Up                1038

num = length(fit.2vs1$genes$Length)
degs.2vs1 <- topTable(fit.2vs1, coef="group2 - group1", confint=TRUE, number = num)
write.table(degs.2vs1, file=paste0('parental','.degs.txt'), sep='\t',row.names = TRUE, quote = FALSE);

############## release VS nilo #################
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
matrix.3vs1 <- makeContrasts(group3-group1,levels=design)
# compute estimated coefficients and standard errors for a given set of contrasts
fit.3vs1 <- contrasts.fit(fit, matrix.3vs1)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit.3vs1 <- eBayes(fit.3vs1)
options(digits=3)

#colnames(fit$coefficients)

summary(decideTests(fit.3vs1, p.value=0.05,lfc=1))
#group3 - group1
#Down              3619
#NotSig           11920
#Up                1226

num = length(fit.3vs1$genes$Length)
degs.3vs1 <- topTable(fit.3vs1, coef="group3 - group1", confint=TRUE, number = num)
write.table(degs.3vs1, file=paste0('release','.degs.txt'), sep='\t',row.names = TRUE, quote = FALSE);

############## release(group3) VS parental(group2) #################NEED RETHINK######!!!!!!!!
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
matrix.3vs2 <- makeContrasts(group3-group2,levels=design)
# compute estimated coefficients and standard errors for a given set of contrasts
fit.3vs2 <- contrasts.fit(fit, matrix.3vs2)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit.3vs2 <- eBayes(fit.3vs2)
options(digits=3)

#colnames(fit$coefficients)

summary(decideTests(fit.3vs2, p.value=0.05,lfc=1))
#group3 - group1
#Down              34
#NotSig           16572
#Up                159

num = length(fit.3vs2$genes$Length)
degs.3vs2 <- topTable(fit.3vs2, coef="group3 - group2", confint=TRUE, number = num)
write.table(degs.3vs2, file=paste0('release-parental','.degs.txt'), sep='\t',row.names = TRUE, quote = FALSE);


###novel lncRNA only###
testgtf <- read.delim("test.combined.gtf", check.names=FALSE,header = FALSE)
testgtf2<-separate(testgtf,V9,c(paste0("V",10:16)),sep=";")
for (i in 1:nrow(testgtf2)){
testgtf2$V17[i]<-strsplit(testgtf2$V11[i]," ")[[1]][3]
}

testgtf3<-testgtf2[testgtf2$V17%in%genelist,]
length(unique(testgtf3$V17)) #219 (novel lncRNA Differential expression)
write.table(testgtf3,"testgtf3.gtf",sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

#####novel lncRNAs##
length(grep("XLOC_",row.names(degs.2vs1),ignore.case=FALSE,fixed=TRUE))
#[1] 219

##novel genes dataset##
novel_2vs1<-degs.2vs1[grep("XLOC_",row.names(degs.2vs1),ignore.case=FALSE,fixed=TRUE),]
novel_3vs1<-degs.3vs1[grep("XLOC_",row.names(degs.3vs1),ignore.case=FALSE,fixed=TRUE),]
novel_3vs2<-degs.3vs2[grep("XLOC_",row.names(degs.3vs2),ignore.case=FALSE,fixed=TRUE),]

novel_2vs1$geneid<-novel_2vs1$genename<-rownames(novel_2vs1)
novel_2vs1<-novel_2vs1[,c(10:11,1:9)]
########known genes (annotated) dataset#########(Jan.2,2019)
known_2vs1<-degs.2vs1[-grep("XLOC_",row.names(degs.2vs1),ignore.case=FALSE,fixed=TRUE),]
known_3vs1<-degs.3vs1[-grep("XLOC_",row.names(degs.3vs1),ignore.case=FALSE,fixed=TRUE),]
known_3vs2<-degs.3vs2[-grep("XLOC_",row.names(degs.3vs2),ignore.case=FALSE,fixed=TRUE),]

row.names(known_2vs1)<-substring(row.names(known_2vs1),1,15)
row.names(known_3vs1)<-substring(row.names(known_3vs1),1,15)
row.names(known_3vs2)<-substring(row.names(known_3vs2),1,15)



########################################################################################################################

#######################################################################################################################
#####2018-08-24 GSEA analysis Data preparation (*.rnk)#################################################################
#######################################################################################################################
#############convert ensemble id to gene symbols#######################
#############install biomaRt package################################
#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
##load biomaRt
library("biomaRt")

ensembl = useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl")
#filters = listFilters(ensembl)
genesym3vs2<-getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), filters = "ensembl_gene_id", values=rownames(known_3vs2), mart=ensembl)

known_3vs2$geneid<-rownames(known_3vs2)
genesym3vs2$geneid<-genesym3vs2$ensembl_gene_id
known_3vs2final<-full_join(known_3vs2,genesym3vs2,by="geneid")
#known_2vs1final$genename<-known_2vs1final$hgnc_symbol
#known_2vs1final<-known_2vs1final[,c(10,13,1:9)]
#allDGE_2vs1<-rbind(known_2vs1final,novel_2vs1)



#######LncRNA list#########
file1<-read.delim("Table_Peaks_in_resistant.gtf",header = FALSE)
file1$V16<-NULL
write.table(file1,"Table_Peaks_in_resistant.gtf",sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

#overlap<-read.delim("testoverlap.txt",header = FALSE) #19 (testgtf3 and Table*.gtf overlap genes)

file2<-read.delim("Table_Peaks_in_Parental.gtf",header = FALSE)
file3<-read.delim("Table_Peaks_in_release.gtf",header = FALSE)

##lncrna files from hg38##
lncrna<-read.delim("lncrna.txt",header=FALSE)
lncrna<-separate(lncrna,V9,c(paste0("V",10:16)),sep=";")
lncrna<-lncrna[lncrna$V3=="transcript",]
for (i in 1:nrow(lncrna)){
    lncrna$V17[i]<-substring(lncrna$V12[i],10,24)
}

for (i in 1:nrow(lncrna)){
    lncrna$V18[i]<-substring(lncrna$V13[i],12)
}

lncrna<-lncrna[lncrna$V16!=" gene_type protein_coding",]


##2vs1 and 1##
known2vs1overlap<-known_2vs1final[known_2vs1final$hgnc_symbol%in%file1$V4,]
known2vs1overlaplncrna<-known2vs1overlap[known2vs1overlap$hgnc_symbol%in%lncrna$V18,] #188
##2vs1 and 2##
known2vs1overlap2<-known_2vs1final[known_2vs1final$hgnc_symbol%in%file2$V4,]
known2vs1overlaplncrna2<-known2vs1overlap2[known2vs1overlap2$hgnc_symbol%in%lncrna$V18,] #162
##1 and 2 diffirent genes##
known2vs1overlaplncrna1_2<-known2vs1overlaplncrna[!(known2vs1overlaplncrna$hgnc_symbol%in%known2vs1overlaplncrna2$hgnc_symbol),]
write.table(known2vs1overlaplncrna1_2,"known2vs1overlaplncrna1_2.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE)



known2vs1overlaplncrna.all<-data.frame(V1=c(setdiff(known2vs1overlaplncrna$hgnc_symbol,known2vs1overlaplncrna2$hgnc_symbol),setdiff(known2vs1overlaplncrna2$hgnc_symbol,known2vs1overlaplncrna$hgnc_symbol)))
known2vs1overlaplncrna.all$V2vs1<-known2vs1overlaplncrna.all$V1
#################

##3vs1 and 1##
known3vs1overlap<-known_3vs1final[known_3vs1final$hgnc_symbol%in%file1$V4,]
known3vs1overlaplncrna<-known3vs1overlap[known3vs1overlap$hgnc_symbol%in%lncrna$V18,] #188
##3vs1 and 2##
known3vs1overlap2<-known_3vs1final[known_3vs1final$hgnc_symbol%in%file3$V4,]
known3vs1overlaplncrna2<-known3vs1overlap2[known3vs1overlap2$hgnc_symbol%in%lncrna$V18,] #162
##1 and 2 diffirent genes##
known3vs1overlaplncrna1_2<-known3vs1overlaplncrna[!(known3vs1overlaplncrna$hgnc_symbol%in%known3vs1overlaplncrna2$hgnc_symbol),]
#known3vs1overlaplncrna2_1<-known3vs1overlaplncrna[!(known3vs1overlaplncrna2$hgnc_symbol%in%known3vs1overlaplncrna$hgnc_symbol),]
write.table(known3vs1overlaplncrna1_2,"known3vs1overlaplncrna1_2.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE)

known3vs1overlaplncrna.all<-data_frame(V1=c(setdiff(known3vs1overlaplncrna$hgnc_symbol,known3vs1overlaplncrna2$hgnc_symbol),setdiff(known3vs1overlaplncrna2$hgnc_symbol,known3vs1overlaplncrna$hgnc_symbol)))
known3vs1overlaplncrna.all$V3vs1<-known3vs1overlaplncrna.all$V1
#################

#################

##3vs2 and 3##
known3vs2overlap<-known_3vs2final[known_3vs2final$hgnc_symbol%in%file3$V4,]
known3vs2overlaplncrna<-known3vs2overlap[known3vs2overlap$hgnc_symbol%in%lncrna$V18,] #188
##3vs2 and 2##
known3vs2overlap2<-known_3vs2final[known_3vs2final$hgnc_symbol%in%file2$V4,]
known3vs2overlaplncrna2<-known3vs2overlap2[known3vs2overlap2$hgnc_symbol%in%lncrna$V18,] #162
##1 and 2 diffirent genes##
known3vs2overlaplncrna1_2<-known3vs2overlaplncrna[!(known3vs2overlaplncrna$hgnc_symbol%in%known3vs2overlaplncrna2$hgnc_symbol),]
write.table(known3vs2overlaplncrna1_2,"known3vs2overlaplncrna1_2.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names = TRUE)


known3vs2overlaplncrna.all<-data_frame(V1=c(setdiff(known3vs2overlaplncrna$hgnc_symbol,known3vs2overlaplncrna2$hgnc_symbol),setdiff(known3vs2overlaplncrna2$hgnc_symbol,known3vs2overlaplncrna$hgnc_symbol)))
known3vs2overlaplncrna.all$V3vs2<-known3vs2overlaplncrna.all$V1
#################


knownalllncrnadiff<-full_join(known2vs1overlaplncrna.all,known3vs1overlaplncrna.all,by="V1")

knownalllncrnadiff<-full_join(knownalllncrnadiff,known3vs2overlaplncrna.all,by="V1")

write.table(knownalllncrnadiff,"knowalllncrnadiff.txt",quote=FALSE,row.names=FALSE,col.names = TRUE,sep="\t")
#########################################################################################################################











