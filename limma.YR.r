#############################Date:2018-08-23#####################################################
#############################RNA-seq analysis-downstream#########################################
#######Input file: counts*.cut.txt files from featureCounts commands in Linux####################
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


# args <- commandArgs(trailingOnly = TRUE);
# input.matrix = args[1];

src.data <- read.delim("countsall.txt", row.names="Geneid", check.names=FALSE);
data <- src.data[,seq(2,10)];
group <- factor(c(1,1,1,2,2,2,3,3,3))


y <- DGEList(counts=data, group=group);
y$genes <- data.frame(Length=src.data$Length);
rownames(y$genes) <- rownames(y$counts);

#y$genes$Symbol <- data.frame(Length=src.data$Symbol);
#rownames(y$symbol) <- rownames(y$counts);

# filtering out low expressed genes
#  the minimum number of samples in each group is three, over here.
keep <- rowSums(cpm(y)>1) >= 3
table(keep);
# FALSE  TRUE 
# 7417 12501
y <- y[keep, keep.lib.sizes=FALSE];



design <- model.matrix(~0+group);
rownames(design) <- colnames(y)
design
#group1 group2 group3
#Vcap_sicrtl_1      1      0      0
#Vcap_sicrtl_2      1      0      0
#Vcap_sicrtl_4      1      0      0
#Vcap_si025_1       0      1      0
#Vcap_si025_2       0      1      0
#Vcap_si025_4       0      1      0
#Vcap_si033_1       0      0      1
#Vcap_si033_2       0      0      1
#Vcap_si033_4       0      0      1
# Note that the filtering does not use knowledge of what treatment corresponds to each sample, so
# the filtering does not bias the subsequent differential expression analysis.
# The TMM normalization is applied to account for the compositional biases:


# normalization by the library sizes
y <- calcNormFactors(y);
y$samples;

write.table(rpkm(y), file=paste0('rpkm','.count.txt'), sep='\t', row.names = TRUE, quote = FALSE);

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
write.table(counts.voom,file="counts.voom.txt",row.names=T,quote=F,sep="\t");
########################################################################################################################
# fit linear model for each gene given a series of libraries
fit <- lmFit(v, design)
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
matrix.2vs1 <- makeContrasts(group2-group1,levels=design)
# compute estimated coefficients and standard errors for a given set of contrasts
fit.2vs1 <- contrasts.fit(fit, matrix.2vs1)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit.2vs1 <- eBayes(fit.2vs1)
options(digits=3)

#colnames(fit$coefficients)

summary(decideTests(fit.2vs1, p.value=0.05,lfc=1))
# group2 - group1
# Down               170
# NotSig           11685
# Up                 646

num = length(fit.2vs1$genes$Length)
degs.2vs1 <- topTable(fit.2vs1, coef="group2 - group1", confint=TRUE, number = num)
write.table(degs.2vs1, file=paste0('si025','.degs.txt'), sep='\t',row.names = TRUE, quote = FALSE);



############## si033 VS control #################
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
matrix.3vs1 <- makeContrasts(group3-group1,levels=design)
# compute estimated coefficients and standard errors for a given set of contrasts
fit.3vs1 <- contrasts.fit(fit, matrix.3vs1)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit.3vs1 <- eBayes(fit.3vs1)
options(digits=3)

#colnames(fit$coefficients)

summary(decideTests(fit.3vs1, p.value=0.05,lfc=1))
# group3 - group1
# Down               241
# NotSig           11912
# Up                 348

num = length(fit.3vs1$genes$Length)
degs.3vs1 <- topTable(fit.3vs1, coef="group3 - group1", confint=TRUE, number = num)
write.table(degs.3vs1, file=paste0('si033','.degs.txt'), sep='\t',row.names = TRUE, quote = FALSE);


#######################################################################################################################
#####2018-08-24 GSEA analysis Data preparation (*.rnk)#################################################################
#######################################################################################################################
results025 <- read.delim("si025.degs.txt", check.names=FALSE)
results033 <- read.delim("si033.degs.txt", check.names=FALSE)
#############convert ensemble id to gene symbols#######################
#############install biomaRt package################################
#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
##load biomaRt
library("biomaRt")
library(dplyr)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)

enidtognf<-function(data){
    test<-rownames(data)
    genesym<-getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), filters = "ensembl_gene_id", values=test, mart=ensembl)
    data<-data[test%in%genesym$ensembl_gene_id,]
    return(data)
}

results033test<-enidtognf(results033)
results033test$geneid<-rownames(results033test)
genesym$geneid<-genesym$ensembl_gene_id
results033final<-full_join(results033test,genesym,by="geneid")
results033final$genename<-results033final$hgnc_symbol


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

#colors <- rep(c("blue","red", "green"), 3)
#plotMDS(y, col=colors[group])


#plotMD(cpm(y, log=TRUE), column=1)
#abline(h=0, col="red", lty=2, lwd=2)

#tmm <- cpm(y$counts);


write.table(rpkm(y), file=paste0('rpkm','.count.txt'), sep='\t', row.names = TRUE, quote = FALSE);
# Inputing RNA-seq counts to clustering or heatmap routines designed for microarray data is not
# straight-forward, and the best way to do this is still a matter of research. To draw a heatmap
# of individual RNA-seq samples, we suggest using moderated log-counts-per-million. This can be
# calculated by cpm with positive values for prior.count, for example

logrpkm <- rpkm(y, prior.count=0, log=TRUE);
write.table(logrpkm, file=paste0('logrpkm','.count.txt'), sep='\t', row.names = TRUE, quote = FALSE);
# where y is the normalized DGEList object. This produces a matrix of log2 counts-per-million
# (logCPM), with undefined values avoided and the poorly defined log-fold-changes for low counts
# shrunk towards zero. Larger values for prior.count produce stronger moderation of the values
# for low counts and more shrinkage of the corresponding log-fold-changes.







fit <- glmQLFit(y, design, robust=TRUE);

# si025 VS control
#qlf.2vs1 <- glmQLFTest(fit, coef=2);
qlf.2vs1 <- glmQLFTest(fit, contrast = c(-1,1,0));
topTags(qlf.2vs1)

# si033 VS control
qlf.3vs1 <- glmQLFTest(fit, contrast = c(-1,0,1));
topTags(qlf.3vs1)

#qlf.32vs1 <- glmQLFTest(fit, contrast = c(-1,0.5,0.5));
#topTags(qlf.32vs1)
print(summary(decideTests(qlf.2vs1,  p.value=0.001)));

print(summary(decideTests(qlf.2vs1,  p.value=0.1, lfc=1)));
print(summary(decideTests(qlf.3vs1,  p.value=0.1, lfc=1)));


#png(filename = paste0(args[1], ".plotMD.png"), width = 2680, height = 1480, units = "px", pointsize = 2, bg = "white",  res = 400)
#plotMD(qlf.2vs1,p.value = 0.1);
#abline(h=c(-1, 1), col="green");
#dev.off();


total <- dim(y$counts)[1]
degs <- topTags(qlf.2vs1, n=total, sort.by="PValue");
#output.colnames <- as.vector(c('Geneid', 'logFC', 'logCPM', 'F', 'PValue', 'FDR'));
#write.table(degs, file='tumor_gtex.degs.txt', sep='\t',row.names = TRUE, quote = FALSE, col.names = output.colnames);
write.table(degs, file=paste0('si025','.degs.txt'), sep='\t',row.names = TRUE, quote = FALSE);


degs <- topTags(qlf.3vs1, n=total, sort.by="PValue");
#output.colnames <- as.vector(c('Geneid', 'logFC', 'logCPM', 'F', 'PValue', 'FDR'));
#write.table(degs, file='tumor_gtex.degs.txt', sep='\t',row.names = TRUE, quote = FALSE, col.names = output.colnames);
write.table(degs, file=paste0('si033','.degs.txt'), sep='\t',row.names = TRUE, quote = FALSE);







#####################################
#####################################
#################Heatmap#########################################
#table<- as.data.frame(degs$table);
#name.list <- rownames(subset(table, FDR<0.01 & (logFC> 2 || logFC < -2)));

#logcpm.matrix <- read.delim("./tumor_gtex.logcpm.txt",row.names="Geneid",check.names=FALSE);

#selected.matrix <- subset(logcpm.matrix, rownames(logcpm.matrix) %in% name.list);

mat <- as.matrix(read.delim("./big_set/LFC.heatmap.txt",row.names="Geneid",check.names=FALSE));
mat <- as.matrix(read.delim("./small_set/LFC.heatmap.txt",row.names="Geneid",check.names=FALSE));
#annotation_col <- data.frame(group=contract$condition)
#rownames(annotation_col) <- contract$Sample_ID;

library("RColorBrewer")
library("pheatmap")

H<-pheatmap(mat,
          show_colnames = T,
          show_rownames = T,
          fontsize_row=5,
          fontsize_col=12,
          cluster_cols = F,
          cluster_rows = T,
          color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(500),
#          #color = colorRampPalette((brewer.pal(n = 11, name ="RdYlBu")))(500),
          #clustering_distance_rows="correlation",
          #clustering_distance_cols="correlation",
#          #annotation_row = annotation_row,
#          annotation_col = annotation_col,
#          #annotation_colors = mat_colors,
           #border_color = NA,
#          scale = "row",
          method="complete")


pheatmap(mat[H$tree_row$order, seq(1,6)],
         show_colnames = T,
         show_rownames = T,
         fontsize_row=3.5,
         fontsize_col=12,
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(500),
         #          #color = colorRampPalette((brewer.pal(n = 11, name ="RdYlBu")))(500),
         #clustering_distance_rows="correlation",
         #clustering_distance_cols="correlation",
         #          #annotation_row = annotation_row,
         #          annotation_col = annotation_col,
         #          #annotation_colors = mat_colors,
         border_color = NA,
         #          scale = "row",
         method="complete")





mat2 <- as.matrix(read.delim("./LFC.degs.lncRNA.heatmap.txt",row.names="Geneid",check.names=FALSE));

#annotation_col <- data.frame(group=contract$condition)
#rownames(annotation_col) <- contract$Sample_ID;


pheatmap(mat2,
         show_colnames = T,
         show_rownames = F,
         fontsize_row=5,
         fontsize_col=12,
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(500),
         #          #color = colorRampPalette((brewer.pal(n = 11, name ="RdYlBu")))(500),
         clustering_distance_rows="correlation",
         clustering_distance_cols="correlation",
         #          #annotation_row = annotation_row,
         #          annotation_col = annotation_col,
         #          #annotation_colors = mat_colors,
         border_color = NA,
         #          scale = "row",
         method="complete")











#################
# dat <- read.table('./known.heatmap.txt', sep="\t", header=TRUE, encoding="UTF-8", row.names = 1);
# mat <- as.matrix(dat);
# png(filename = "known.heatmap.png",
#     width = 2680, height = 1480, units = "px", pointsize = 2,
#     bg = "white",  res = 400)
# 
# H <- pheatmap(mat,
#               show_colnames = T,
#               show_rownames = F,
#               fontsize_row=6,
#               fontsize_col=6,
#               cluster_cols = T,
#               cluster_rows = T,
#               color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(500),
#               #color = colorRampPalette((brewer.pal(n = 11, name ="RdYlBu")))(500),
#               clustering_distance_rows="correlation",
#               #clustering_distance_cols="correlation",
#               #annotation_row = annotation_row,
#               #annotation_col = annotation_col,
#               annotation_colors = mat_colors,
#               border_color = NA,
#               scale = "row",
#               method="complete")
# 
# 
# 
# pheatmap(mat[H$tree_row$order, H$tree_col$order],
#          show_colnames = T,
#          fontsize_row=0.1,
#          fontsize_col=2,
#          cluster_rows = F,
#          color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(500),
#          clustering_distance_rows="correlation",
#          clustering_distance_cols="correlation",
#          #annotation_row = annotation_row,
#          #annotation_col = annotation_col,
#          border_color = NA,
#          scale = "row",
#          method="complete")
# dev.off()
# 
