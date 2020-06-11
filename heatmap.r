##########################################################################
#################Heatmap for yuanguo wang 2019-02-04######################
##########################################################################
#table<- as.data.frame(degs$table);
#name.list <- rownames(subset(table, FDR<0.01 & (logFC> 2 || logFC < -2)));
#logcpm.matrix <- read.delim("./tumor_gtex.logcpm.txt",row.names="Geneid",check.names=FALSE);
#selected.matrix <- subset(logcpm.matrix, rownames(logcpm.matrix) %in% name.list);
##yuanguo data##
dat<-read.delim("ygwang.csv",header = FALSE,sep=",")
row.names(dat)<-dat$V1
dat$V1<-NULL
#calculate Z-score by each gene(row)##
zscore<-matrix(data = NA,nrow=5,ncol = 6)
for (i in 1:nrow(dat))
{
    zscore[i,1]<-(dat[i,1]-rowMeans(dat)[i])/sd(dat[i,])
    zscore[i,2]<-(dat[i,2]-rowMeans(dat)[i])/sd(dat[i,])
    zscore[i,3]<-(dat[i,3]-rowMeans(dat)[i])/sd(dat[i,])
    zscore[i,4]<-(dat[i,4]-rowMeans(dat)[i])/sd(dat[i,])
    zscore[i,5]<-(dat[i,5]-rowMeans(dat)[i])/sd(dat[i,])
    zscore[i,6]<-(dat[i,6]-rowMeans(dat)[i])/sd(dat[i,])
}
row.names(zscore)<-row.names(dat)
colnames(zscore)<-paste0(rep(c("V","T"),each=3),1:3)
#calculte z-score by all genes and samples##
avg<-sum(dat)/(5*6) #563.2106
vector<-c(dat$V2,dat$V3,dat$V4,dat$V5,dat$V6,dat$V7)
std<-sd(vector) #1081.213
zscore<-matrix(data = NA,nrow=5,ncol = 6)
for (i in 1:nrow(dat))
{
    zscore[i,1]<-(dat[i,1]-avg)/std
    zscore[i,2]<-(dat[i,2]-avg)/std
    zscore[i,3]<-(dat[i,3]-avg)/std
    zscore[i,4]<-(dat[i,4]-avg)/std
    zscore[i,5]<-(dat[i,5]-avg)/std
    zscore[i,6]<-(dat[i,6]-avg)/std
}
row.names(zscore)<-row.names(dat)
colnames(zscore)<-paste0(rep(c("V","T"),each=3),1:3)

#calculate z-score by samples(column)#
zscore<-matrix(data = NA,nrow=5,ncol = 6)
for (i in 1:ncol(dat))
{
    zscore[1,i]<-(dat[1,i]-colMeans(dat)[i])/sd(dat[,i])
    zscore[2,i]<-(dat[2,i]-colMeans(dat)[i])/sd(dat[,i])
    zscore[3,i]<-(dat[3,i]-colMeans(dat)[i])/sd(dat[,i])
    zscore[4,i]<-(dat[4,i]-colMeans(dat)[i])/sd(dat[,i])
    zscore[5,i]<-(dat[5,i]-colMeans(dat)[i])/sd(dat[,i])
    #zscore[6,i]<-(dat[6,i]-colMeans(dat)[i])/sd(dat[,i])
}
row.names(zscore)<-row.names(dat)
colnames(zscore)<-paste0(rep(c("V","T"),each=3),1:3)


##LUKE data to test if zscore are the same##
dat<-read.delim("gene_heatmap.csv",header = TRUE,sep=",")
row.names(dat)<-dat$Stiffness.kPa.
dat$Stiffness.kPa.<-NULL
#calculate Z-score by each gene##
zscore<-matrix(data = NA,nrow=5,ncol = 5)
for (i in 1:nrow(dat))
{
    zscore[i,1]<-(dat[i,1]-rowMeans(dat)[i])/sd(dat[i,])
    zscore[i,2]<-(dat[i,2]-rowMeans(dat)[i])/sd(dat[i,])
    zscore[i,3]<-(dat[i,3]-rowMeans(dat)[i])/sd(dat[i,])
    zscore[i,4]<-(dat[i,4]-rowMeans(dat)[i])/sd(dat[i,])
    zscore[i,5]<-(dat[i,5]-rowMeans(dat)[i])/sd(dat[i,])
    #zscore[i,6]<-(dat[i,6]-rowMeans(dat)[i])/sd(dat[i,])
}


#mat<-as.matrix(dat)


#mat <- as.matrix(read.delim("./big_set/LFC.heatmap.txt",row.names="Geneid",check.names=FALSE));
#mat <- as.matrix(read.delim("./small_set/LFC.heatmap.txt",row.names="Geneid",check.names=FALSE));
#annotation_col <- data.frame(group=contract$condition)
#rownames(annotation_col) <- contract$Sample_ID;

library("RColorBrewer")
library("pheatmap")

#png(filename = "yuanguo.heatmap.png",width = 2680, height = 1480, units = "px", pointsize = 2,bg = "white",  res = 400)
tiff(filename = "yuanguo.heatmap_zscore(sample).tiff",width = 4.5, height = 3,units='in',res = 400,pointsize = 4)
H<-pheatmap(zscore,
            show_colnames = T,
            show_rownames = T,
            fontsize_row=12,
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
            #border_color = NA,
            #          scale = "row",
            #cellwidth = 15, 
            #cellheight = 12,
            method="complete")
dev.off()



##darker color for ygwang##
png(filename = "yuanguo.heatmap_darker_LUKEtest.png",width = 4, height = 3,units='in',res = 400,pointsize = 4)
pheatmap(zscore,
         show_colnames = T,
         show_rownames = T,
         fontsize_row=6,
         fontsize_col=6,
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(rev(brewer.pal(n =11, name ="RdYlBu")))(500),
         #color = colorRampPalette((brewer.pal(n = 11, name ="RdYlBu")))(500),
         #clustering_distance_rows="euclidean",
         clustering_distance_cols="correlation",
         #annotation_row = annotation_row,
         #annotation_col = annotation_col,
         #annotation_colors = mat_colors,
         #border_color = NA,
         #scale = "row",
         treeheight_row = 0,
         method="complete"
         #main = "log RPKM for DEGs"
         )
dev.off()







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
