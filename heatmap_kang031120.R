###########2020-03-11##################
##heatmap zscore by gene and samples###
#######################################

dat<-read.delim("Kang_03112020.csv",header = TRUE,sep=",",stringsAsFactors = FALSE,check.names = FALSE)

# log2 #
dat[dat == 0] <- 1
dat<-log2(dat)

# gene & sample zscore # both_zscore function from heatmap_function.R #
zscore_both<-both_zscore(dat)

heatmapR("Dr.kang.log2.zscoreboth.tiff",zscore_both,8,8)