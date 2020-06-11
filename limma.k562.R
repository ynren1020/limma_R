#############################Date:2018-11-27#####################################################
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


# args <- commandArgs(trailingOnly = TRUE);
# input.matrix = args[1];

nilo <- read.delim("nilo_id_counts.txt", row.names="Geneid", check.names=FALSE)
nilosub<-nilo[,c(5,6,7)]
nilosub<-rename(nilosub,nilo1="/home/ryang/project/k562_drug/nilo1_IN_sorted.bam")
nilosub<-rename(nilosub,nilo2="/home/ryang/project/k562_drug/nilo2_IN_sorted.bam")

##FUNCTION TO READ SAMPLE DATA IN##
sample<-function(file){
    temp<-read.delim(file,row.names="Geneid",check.names=FALSE)
    tempsub<-temp[,c(5,6,7)]
    colnames(tempsub)[2]<-strsplit(strsplit(colnames(tempsub)[2],"/")[[1]][6],"_")[[1]][1]
    colnames(tempsub)[3]<-strsplit(strsplit(colnames(tempsub)[3],"/")[[1]][6],"_")[[1]][1]
    return(tempsub)
    }
nilo<-sample("nilo_id_counts.txt")
parental<-sample("parental_id_counts.txt")
release<-sample("release_id_counts.txt")

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
#44337 16907 
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
#Down              3347
#NotSig           12511
#Up                1049

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
#Down              3708
#NotSig           11956
#Up                1243

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
#NotSig           16711
#Up                162

num = length(fit.3vs2$genes$Length)
degs.3vs2 <- topTable(fit.3vs2, coef="group3 - group2", confint=TRUE, number = num)
write.table(degs.3vs2, file=paste0('release-parental','.degs.txt'), sep='\t',row.names = TRUE, quote = FALSE);

#####novel genes##
length(grep("G.",row.names(degs.2vs1),ignore.case=FALSE,fixed=TRUE))
#[1] 402

##novel genes dataset##
novel_2vs1<-degs.2vs1[grep("G.",row.names(degs.2vs1),ignore.case=FALSE,fixed=TRUE),]
novel_3vs1<-degs.3vs1[grep("G.",row.names(degs.3vs1),ignore.case=FALSE,fixed=TRUE),]
novel_3vs2<-degs.3vs2[grep("G.",row.names(degs.3vs2),ignore.case=FALSE,fixed=TRUE),]

novel_2vs1$geneid<-novel_2vs1$genename<-rownames(novel_2vs1)
novel_2vs1<-novel_2vs1[,c(10:11,1:9)]
##known genes dataset##
known_2vs1<-degs.2vs1[-grep("G.",row.names(degs.2vs1),ignore.case=FALSE,fixed=TRUE),]
known_3vs1<-degs.3vs1[!grep("G.",row.names(degs.3vs1),ignore.case=FALSE,fixed=TRUE),]
known_3vs2<-degs.3vs2[!grep("G.",row.names(degs.3vs2),ignore.case=FALSE,fixed=TRUE),]

row.names(known_2vs1)<-substring(row.names(known_2vs1),1,15)


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
library(dplyr)
ensembl = useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl")
#filters = listFilters(ensembl)
genesym2vs1<-getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), filters = "ensembl_gene_id", values=rownames(known_2vs1), mart=ensembl)

known_2vs1$geneid<-rownames(known_2vs1)
genesym2vs1$geneid<-genesym2vs1$ensembl_gene_id
known_2vs1final<-full_join(known_2vs1,genesym2vs1,by="geneid")
known_2vs1final$genename<-known_2vs1final$hgnc_symbol
known_2vs1final<-known_2vs1final[,c(10,13,1:9)]


allDGE_2vs1<-rbind(known_2vs1final,novel_2vs1)




#########################################################################################################################

#Create ranks:
rankf<-function(results){
    results.ord <- results[ order(-results[,"logFC"]), ]
    results.ord<-filter(results.ord,genename!="")
    ranks<-data.frame("genename"=results.ord[,ncol(results.ord)],"rank"=results.ord[,"logFC"])
    #ranks<-data.frame("rank"=results.ord[,"logFC"])
    #names(ranks) <- results.ord[,ncol(results.ord)]
   return(ranks)
    
}

rank033<-rankf(results=results033final)
#write.table(rank033, file=paste0('rank033','.txt.rnk'), sep='\t',row.names = TRUE, quote = FALSE);
#write.table(rank033, file=paste0('rank033','.csv.rnk'), sep=',',row.names = TRUE, quote = FALSE);
write.table(rank033, file=paste0('rank033','.txt'), sep='\t',row.names = FALSE, col.names=FALSE,quote = FALSE);



#####################################################################################################################################
#####################################################################################################################################

