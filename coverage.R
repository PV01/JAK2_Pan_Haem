#################
# JAK2 V617F query ##

setwd("H:/R_projects/JAK2_Pan_Haem")

library(tidyverse)
library(odbc)
library(Rsamtools)
library(GenomicRanges)


NGS_coverage<-function(bam_file, chromosome, start, end){
  
  nuc_flag<-ifelse(start==end, TRUE, FALSE)
  
  rl <- IRangesList()
  rl[[chromosome]] <- IRanges(start=start, end=end, names=chromosome)
  sbp <- ScanBamParam(which=rl)
  pp <- PileupParam(max_depth=10^6,
                    min_base_quality=1,
                    min_mapq=5,
                    min_nucleotide_depth=1,
                    distinguish_strands=FALSE,
                    distinguish_nucleotides=nuc_flag,
                    ignore_query_Ns=FALSE)
  p <- pileup(file=bam_file, yieldSize=10^6, scanBamParam=sbp, pileupParam=pp)
  
  
  if(start!=end){
    p<-p %>%
        mutate(mean_coverage=round(mean(count),0),
               pos=which_label) %>%
        distinct(seqnames, pos, mean_coverage)
  }
  
  return(p)
  
}

################################################################


NGS_path<-"S:/MiSeq_data/GRCh38/TSMP_Flex/pipeline_output/"
NGS_path<-list.files(NGS_path, full.names=TRUE, recursive=TRUE)
NGS_path<-NGS_path[!grepl("pre-go-live", NGS_path)]

bam_files<-NGS_path[grepl("sorted.bam$", NGS_path)]
bam_files<-bam_files[!grepl("NEG_S1", bam_files)]

saveRDS(bam_files,"bam_files.RDS")



coverage<-data.frame()
chromosome<-"chr9"
start<-5073770
end<-5073770

for (i in 1: length(bam_files)){
  bam_coverage<-NGS_coverage(bam_files[i], chromosome, start, end)%>%
                            mutate(Worksheet=as.numeric(gsub("TSMP_","",str_extract(bam_files[i],"TSMP_[0-9]{6}"))),
                                  LABNO = gsub("-",".",str_extract(bam_files[i],"D[0-9]{2}-[0-9]{5}")))
  coverage<-rbind(coverage,bam_coverage)
  
  print(i)
              
}

saveRDS(coverage,"coverage.RDS")
