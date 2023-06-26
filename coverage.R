#################
# JAK2 V617F query ##

setwd("H:/R_projects/JAK2_Pan_Haem")

library(tidyverse)
library(odbc)
library(Rsamtools)
library(GenomicRanges)
library(ggplot2)


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

#saveRDS(coverage,"coverage.RDS")


coverage<-readRDS("coverage.RDS")


coverage_summary<-coverage %>%
                          mutate(nuc_count=paste0(nucleotide,":",count)) %>%
                          group_by(Worksheet, LABNO) %>%
                          mutate(calls=paste0(nuc_count, collapse="/"),
                                 Total=sum(count)) %>%
                          distinct(Worksheet, LABNO, .keep_all = TRUE) %>%
                          select(LABNO, Worksheet, calls, Total)

mean_depth<-mean(coverage_summary$Total)

ggplot(data=coverage_summary, aes(x="A",y=Total))+
      geom_boxplot()+
      geom_jitter(alpha = 0.3, color = "darkblue")


m<-mean(coverage_summary$Total)
s<-sd(coverage_summary$Total)


coverage_summary_Z<-coverage_summary %>%
                       mutate(Z=(Total-m)/s) %>%
                       filter(Z>-3, Z<3)
                
summary(10/coverage_summary_Z$Total*100)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.793    1.319    1.610    2.423    2.288 1000.000 
hist(10/coverage_summary_Z$Total*1000)




ggplot(data=coverage_summary, aes(x="A",y=(10/Total)*100))+
  geom_boxplot(outlier.shape = NA)+
  ylim(0,5)




