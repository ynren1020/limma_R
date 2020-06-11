##2019-02-13##
##create splitted pbs file##
############################

splitfile<-c("xab","xac","xad","xae","xaf","xag","xah")
for (i in 1:length(splitfile))
{
    sink(paste0(splitfile[i],'.pbs'))
    cat("#!/bin/bash -l \n")
    cat("#PBS -l walltime=12:00:00,nodes=1:ppn=8,pmem=1000mb \n")
    cat("#PBS -m abe \n")
    cat("#PBS -M renxx275@umn.edu \n")
    cat("cd $PBS_O_WORKDIR \n")
    cat("module load hisat2 \n")
    cat("module load samtools\n")
    cat(paste0("bash ",splitfile[i], "\n"))
    sink()
}

