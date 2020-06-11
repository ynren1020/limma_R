##########################################################################
#################Heatmap Function 2019-06-20######################
##########################################################################
##load packages##
library("RColorBrewer")
library("pheatmap")
library(colorspace)
library(colorRamps)
library(grDevices)

##Function1.load df and make matrix and take log10 transformation##
create_matrix<-function(input,header){
dat<-read.delim(input,header = header,sep=",",stringsAsFactors = FALSE)
row.names(dat)<-dat$ID
dat$ID<-NULL
dat[dat == 0] <- 1
dat<-log2(dat)
return(dat)
}

dat<-create_matrix("kang.csv",TRUE)

##function2. zscore by genes##
gene_zscore<-function(dat){
#calculate Z-score by each gene(row)##
zscore<-matrix(data = NA,nrow=dim(dat)[1],ncol = dim(dat)[2])
for (i in 1:nrow(dat))
{for (j in 1:ncol(dat)){
    zscore[i,j]<-(dat[i,j]-rowMeans(dat)[i])/sd(dat[i,])
    
}
}
row.names(zscore)<-row.names(dat)
colnames(zscore)<-colnames(dat)
return(zscore)
}

##apply gene_score function##
zscore_gene<-gene_zscore(dat)

##Function3.zscore calculation by sample and genes, 2 dimensions##
both_zscore<-function(dat){
#calculte z-score by all genes and samples## USE THIS ONE!
#AFTER LOG TRANSFORM, THE VARIATION WAS GREATLY REDUCED AND HEATMAP LOOKS BETTER##
avg<-sum(dat)/(dim(dat)[1]*dim(dat)[2]) #3.2459
#vector<-c(dat$DMSO.1,dat$DMSO.2,dat$DMSO.3,dat$DMSO.TGF..1,dat$DMSO.TGF..2,dat$DMSO.TGF..3,dat$PF228.1,dat$PF228.2,dat$PF228.3,dat$PF228.TGF..1,dat$PF228.TGF..2,dat$PF228.TGF..3)
#df to vector by column
vector<-unlist(dat)
std<-sd(vector) # 1.223694
zscore<-matrix(data = NA,nrow=dim(dat)[1],ncol = dim(dat)[2])
for (i in 1:nrow(dat))
{for (j in 1:ncol(dat)){
    zscore[i,j]<-(dat[i,j]-avg)/std
    
}}
row.names(zscore)<-row.names(dat)
colnames(zscore)<-colnames(dat)
return(zscore)
}
##apply both_zscore to normalize expression by all values
zscore_both<-both_zscore(dat)
######################################
##Function4.zscore by samples#########
sample_zscore<-function(dat){
#calculate z-score by samples(column)#
zscore<-matrix(data = NA,nrow=dim(dat)[1],ncol=dim(dat)[2])
for (i in 1:ncol(dat))
 { 
    for (j in 1:nrow(dat)){
    zscore[j,i]<-(dat[j,i]-colMeans(dat)[i])/sd(dat[,i])
    
    }
}
row.names(zscore)<-row.names(dat)
colnames(zscore)<-colnames(dat)
return(zscore)
}

zscore_sample<-sample_zscore(dat)



heatmapR<-function(name,zscore,sizeRow,sizeCol){
#png(filename = "yuanguo.heatmap.png",width = 2680, height = 1480, units = "px", pointsize = 2,bg = "white",  res = 400)
tiff(filename = name,width = 4.5, height = 3,units='in',res = 400,pointsize = 4)
H<-pheatmap(zscore,
            show_colnames = T,
            show_rownames = F,
            #breaks = myBreaks,
            fontsize_row=sizeRow,
            fontsize_col=sizeCol,
            cluster_cols = F,
            cluster_rows = F,
            #color = diverge_hcl(12,h=c(246,40)),
            color = matlab.like2(500),
            #color = terrain.colors(500),
            #color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(10000-1),
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

}
heatmapR("Dr.kang.log2.gene.matlab.tiff",zscore_gene,8,8)

paletteLength<-10000
myBreaks <- c(seq(min(zscore_gene), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(zscore_gene)/paletteLength, max(zscore_gene), length.out=floor(paletteLength/2)))



myColor <- colorRampPalette(c("Brown","yellow", "white", "deepskyblue","blue","purple"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths

myColor=matlab.like2(500)
# Plot the heatmap
pheatmap(zscore_gene, color=myColor)
