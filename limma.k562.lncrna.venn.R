###############2020-06-11################
##lncRNAs' venn diagram for Dr. Liu #####
#########################################

degs.2vs1 <- data.table::fread("parental.degs.txt") #16765

known_2vs1<-degs.2vs1[-grep("XLOC_",degs.2vs1$V1,ignore.case=FALSE,fixed=TRUE),] #16546

known_2vs1$gene_id<-substring(known_2vs1$V1,1,15)

# gene id to gene symbol ---

library("biomaRt")

ensembl = useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl")
#filters = listFilters(ensembl)
genesym2vs1<-getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), filters = "ensembl_gene_id", values=known_2vs1$gene_id, mart=ensembl)

genesym2vs1$gene_id<-genesym2vs1$ensembl_gene_id
known_2vs1final<-dplyr::full_join(known_2vs1,genesym2vs1,by="gene_id") #16546
known_2vs1final$genename<-known_2vs1final$hgnc_symbol
known_2vs1final<-known_2vs1final[,c(11,14,2:10)]

# m6a ---
file1 <- data.table::fread("Table_Peaks_in_resistant.gtf")
file2 <- data.table::fread("Table_Peaks_in_Parental.gtf")

##lncrna files from hg38##
lncrna<-read.delim("lncrna.txt",header=FALSE)
lncrna<-tidyr::separate(lncrna,V9,c(paste0("V",10:16)),sep=";")
lncrna<-lncrna[lncrna$V3=="transcript",]

lncrna<-read.delim("lncrna.txt",header=FALSE)
lncrna<-tidyr::separate(lncrna,V9,c(paste0("V",10:16)),sep=";")
lncrna<-lncrna[lncrna$V3=="transcript",] #52774

lncrna<-lncrna[lncrna$V16!=" gene_type protein_coding" & lncrna$V15!=" gene_type protein_coding",] #51586

any(grepl("gene_name", lncrna$V12)) #FALSE
any(grepl("gene_name", lncrna$V13)) #TRUE
any(grepl("gene_name", lncrna$V14)) #FALSE

for (i in 1:nrow(lncrna)){
    lncrna$V17[i]<-substring(lncrna$V12[i],10,24)
    lncrna$V18[i]<-substring(lncrna$V13[i],12)
    #lncrna$V18[i]<-strsplit(lncrna$V13[i], " ")[[1]][3]
}



##how many known lncRNA differential expressed##Jan 7th 2019##
known2vs1lncrna<-known_2vs1final[known_2vs1final$gene_id%in%lncrna$V17,] # 2213
known2vs1lncrna<-known_2vs1final[known_2vs1final$genename%in%lncrna$V18,] #539

write.table(known2vs1lncrna, "known2vs1lncrna_update.txt", sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
##nilo high
known2vs1lncrna_negative<-known2vs1lncrna[known2vs1lncrna$logFC<0&known2vs1lncrna$adj.P.Val<0.05,] 
known2vs1lncrna_negative <- known2vs1lncrna_negative[known2vs1lncrna_negative$genename!=""&(!is.na(known2vs1lncrna_negative$genename)),] ##315 -->nilo high expressed
known2vs1lncrna_negative_1<-known2vs1lncrna_negative[known2vs1lncrna_negative$genename%in%file1$V4&!(known2vs1lncrna_negative$genename%in%file2$V4),] #40 #only in resistant not in parental
known2vs1lncrna_negative_1<-dplyr::arrange(known2vs1lncrna_negative_1, logFC)
write.table(known2vs1lncrna_negative_1, "known2vs1lncrna_negative_1_update.txt", sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
write.table(known2vs1lncrna_negative, "known2vs1lncrna_negative_update.txt", sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)


##parental high
known2vs1lncrna_positive<-known2vs1lncrna[known2vs1lncrna$logFC>0&known2vs1lncrna$adj.P.Val<0.05,] 
known2vs1lncrna_positive <- known2vs1lncrna_positive[known2vs1lncrna_positive$genename!=""&(!is.na(known2vs1lncrna_positive$genename)),] #97 -->parental high expressed
known2vs1lncrna_positive_1<-known2vs1lncrna_positive[known2vs1lncrna_positive$genename%in%file2$V4&!(known2vs1lncrna_positive$genename%in%file1$V4),] #12 #only in parental not in resistant
known2vs1lncrna_positive_1 <- dplyr::arrange(known2vs1lncrna_positive_1, desc(logFC))
write.table(known2vs1lncrna_positive_1, "known2vs1lncrna_positive_1_update.txt", sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
write.table(known2vs1lncrna_positive, "known2vs1lncrna_positive_update.txt", sep="\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)


# output to xlsx ---
library(xlsx)
write.xlsx(known2vs1lncrna_negative_1, file="lncRNAs_m6a.xlsx", sheetName="DElncRNAs_resistance_m6a", row.names=FALSE)
write.xlsx(known2vs1lncrna_negative, file="lncRNAs_m6a.xlsx", sheetName="DElncRNAs_resistance", append=TRUE, row.names=FALSE)
write.xlsx(known2vs1lncrna_positive_1, file="lncRNAs_m6a.xlsx", sheetName="DElncRNAs_parental_m6a", append=TRUE,row.names=FALSE)
write.xlsx(known2vs1lncrna_positive, file="lncRNAs_m6a.xlsx", sheetName="DElncRNAs_parental", append=TRUE, row.names=FALSE)


# Venn ---
library(VennDiagram)

# nilo high ---
set1 <- known2vs1lncrna_negative$genename
set2 <- unique(file1$V4)
set3 <- unique(file2$V4)

# parental high ---
set1 <- known2vs1lncrna_positive$genename
set2 <- unique(file1$V4)
set3 <- unique(file2$V4)



# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
    x = list(set1, set2, set3),
    category.names = c("Highly expressed \n in parental" , "m6a resistance " , "m6a parental"),
    filename = 'lncrna_highinParental_venn_diagramm.png',
    output=TRUE,
    
    # Output features
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.3,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
)

# heatmap ---
library("RColorBrewer")
library("pheatmap")

lncrnaneg1 <- known2vs1lncrna_negative_1
#lncrnaneg1<-read.delim("known2vs1lncrna_negative_1_update.txt",sep="\t")
lncrnaneg1sub<-lncrnaneg1[,c("logFC","genename")]
.rowNamesDF(lncrnaneg1sub,make.names=TRUE)<-lncrnaneg1sub$genename
lncrnaneg1sub$logFC<-abs(lncrnaneg1sub$logFC)
#lncrnaneg1sub<-lncrnaneg1sub[order(-lncrnaneg1sub$logFC),]
lncrnaneg1sub$genename<-NULL
lncrnaneg1sub_mat<-as.matrix(lncrnaneg1sub)
rownames(lncrnaneg1sub_mat) <- lncrnaneg1$genename
#lncrnaneg1suborder<-lncrnaneg1sub[order(lncrnaneg1sub$logFC),]


png(filename = "logFC.heatmap.red6_update.png",width = 3, height = 6,units='in',res = 400,pointsize = 4)

pheatmap(lncrnaneg1sub_mat,
         show_colnames = F,
         show_rownames = T,
         fontsize_row=10,
         fontsize_col=12,
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(brewer.pal(n = 11, name ="Reds"))(500),
         #          #color = colorRampPalette((brewer.pal(n = 11, name ="RDYlBu")))(500),
         #clustering_distance_rows="correlation",
         #clustering_distance_cols="correlation",
         #          #annotation_row = annotation_row,
         #          annotation_col = annotation_col,
         #          #annotation_colors = mat_colors,
         border_color = NA,
         #          scale = "row",
         width = 1,
         height = 4,
         method="complete")

dev.off()

