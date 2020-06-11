############Nov 26th 2018###########
######merge cp_nilo.gtf,parental.annotated.gtf,release.annotated.gtf######
######ref : cp_nilo.gtf######

##load packages##
library(dplyr)
library(tidyr)
library(stringr)


##read data##
ref<-read.delim("cp_nilo.gtf",header = FALSE)
p.anno<-read.delim("parental.annotated.gtf",header = FALSE)
r.anno<-read.delim("release.annotated.gtf",header = FALSE)

##reshape data##
p.anno1<-separate(p.anno,V9,c("t_id","g_id","g_name","xloc","cmp_ref","class_code","tss_id"),sep = ";")
p.anno2<-filter(p.anno1,class_code==" class_code =")

##1 for loop for change ids##works#
for (i in 1:nrow(p.anno1)){
    
if (p.anno1$V3[i]=="transcript"){
    if (p.anno1$t_id[i]%in%p.anno2$t_id) {p.anno1$t_id_new[i]<-str_replace(p.anno1$t_id[i],str_split(p.anno1$t_id[i]," ")[[1]][2],str_split(p.anno1$g_name[i]," ")[[1]][3])
                                          p.anno1$g_id_new[i]<-str_replace(p.anno1$g_id[i],str_split(p.anno1$g_id[i]," ")[[1]][3],str_split(p.anno1$g_name[i]," ")[[1]][3]) }
    else{p.anno1$t_id_new[i]<-paste0(p.anno1$t_id[i],".p")
         p.anno1$g_id_new[i]<-paste0(p.anno1$g_id[i],".p")}

}else if (p.anno1$V3[i]=="exon"){
    
    if (p.anno1$t_id[i]%in%p.anno2$t_id) {p.anno1$t_id_new[i]<-str_replace(p.anno1$t_id[i],str_split(p.anno1$t_id[i]," ")[[1]][2],str_split(p.anno1$t_id_new[i-1]," ")[[1]][2])
                                          p.anno1$g_id_new[i]<-str_replace(p.anno1$g_id[i],str_split(p.anno1$g_id[i]," ")[[1]][3],str_split(p.anno1$g_id_new[i-1]," ")[[1]][3]) }
    else{p.anno1$t_id_new[i]<-paste0(p.anno1$t_id[i],".p")
         p.anno1$g_id_new[i]<-paste0(p.anno1$g_id[i],".p")}
        
} else {print ("error!")}

    
}

##2 for loop for change ids##works#
for (i in 1:nrow(p.anno1)){
    
    
        if (p.anno1$t_id[i]%in%p.anno2$t_id & p.anno1$V3[i]=="transcript") {p.anno1$t_id_new[i]<-str_replace(p.anno1$t_id[i],str_split(p.anno1$t_id[i]," ")[[1]][2],str_split(p.anno1$g_name[i]," ")[[1]][3])
        p.anno1$g_id_new[i]<-str_replace(p.anno1$g_id[i],str_split(p.anno1$g_id[i]," ")[[1]][3],str_split(p.anno1$g_name[i]," ")[[1]][3]) 
        }else if (p.anno1$t_id[i]%in%p.anno2$t_id & p.anno1$V3[i]=="exon"){p.anno1$t_id_new[i]<-str_replace(p.anno1$t_id[i],str_split(p.anno1$t_id[i]," ")[[1]][2],str_split(p.anno1$t_id_new[i-1]," ")[[1]][2])
        p.anno1$g_id_new[i]<-str_replace(p.anno1$g_id[i],str_split(p.anno1$g_id[i]," ")[[1]][3],str_split(p.anno1$g_id_new[i-1]," ")[[1]][3])
        }else{p.anno1$t_id_new[i]<-paste0(p.anno1$t_id[i],".p")
              p.anno1$g_id_new[i]<-paste0(p.anno1$g_id[i],".p")}
    
}

##remove columns##
p.anno3<-p.anno1[,c(1:8,17,16,11)]
p.anno3<-unite(p.anno3,V9,c("g_id_new","t_id_new","g_name"),sep = ";")
write.table(p.anno3,file="cp_parental_new.gtf",row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t")

###############release gtf#####################
##reshape data##
r.anno1<-separate(r.anno,V9,c("t_id","g_id","g_name","xloc","cmp_ref","class_code","tss_id"),sep = ";")
r.anno2<-filter(r.anno1,class_code==" class_code =")

##1 for loop for change ids##works#
for (i in 1:nrow(r.anno1)){
    
    if (r.anno1$V3[i]=="transcript"){
        if (r.anno1$t_id[i]%in%r.anno2$t_id) {r.anno1$t_id_new[i]<-str_replace(r.anno1$t_id[i],str_split(r.anno1$t_id[i]," ")[[1]][2],str_split(r.anno1$g_name[i]," ")[[1]][3])
        r.anno1$g_id_new[i]<-str_replace(r.anno1$g_id[i],str_split(r.anno1$g_id[i]," ")[[1]][3],str_split(r.anno1$g_name[i]," ")[[1]][3]) }
        else{r.anno1$t_id_new[i]<-paste0(r.anno1$t_id[i],".r")
        r.anno1$g_id_new[i]<-paste0(r.anno1$g_id[i],".r")}
        
    }else if (r.anno1$V3[i]=="exon"){
        
        if (r.anno1$t_id[i]%in%r.anno2$t_id) {r.anno1$t_id_new[i]<-str_replace(r.anno1$t_id[i],str_split(r.anno1$t_id[i]," ")[[1]][2],str_split(r.anno1$t_id_new[i-1]," ")[[1]][2])
        r.anno1$g_id_new[i]<-str_replace(r.anno1$g_id[i],str_split(r.anno1$g_id[i]," ")[[1]][3],str_split(r.anno1$g_id_new[i-1]," ")[[1]][3]) }
        else{r.anno1$t_id_new[i]<-paste0(r.anno1$t_id[i],".r")
        r.anno1$g_id_new[i]<-paste0(r.anno1$g_id[i],".r")}
        
    } else {print ("error!")}
    
    
}

##2 for loop for change ids##works#
for (i in 1:nrow(r.anno1)){
    
    
    if (r.anno1$t_id[i]%in%r.anno2$t_id & r.anno1$V3[i]=="transcript") {r.anno1$t_id_new[i]<-str_replace(r.anno1$t_id[i],str_split(r.anno1$t_id[i]," ")[[1]][2],str_split(r.anno1$g_name[i]," ")[[1]][3])
    r.anno1$g_id_new[i]<-str_replace(r.anno1$g_id[i],str_split(r.anno1$g_id[i]," ")[[1]][3],str_split(r.anno1$g_name[i]," ")[[1]][3]) 
    }else if (r.anno1$t_id[i]%in%r.anno2$t_id & r.anno1$V3[i]=="exon"){r.anno1$t_id_new[i]<-str_replace(r.anno1$t_id[i],str_split(r.anno1$t_id[i]," ")[[1]][2],str_split(r.anno1$t_id_new[i-1]," ")[[1]][2])
    r.anno1$g_id_new[i]<-str_replace(r.anno1$g_id[i],str_split(r.anno1$g_id[i]," ")[[1]][3],str_split(r.anno1$g_id_new[i-1]," ")[[1]][3])
    }else{r.anno1$t_id_new[i]<-paste0(r.anno1$t_id[i],".r")
    r.anno1$g_id_new[i]<-paste0(r.anno1$g_id[i],".r")}
    
}

##remove columns##
r.anno3<-r.anno1[,c(1:8,17,16,11)]
r.anno3<-unite(r.anno3,V9,c("g_id_new","t_id_new","g_name"),sep = ";")
write.table(r.anno3,file="cp_release_new.gtf",row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t")




