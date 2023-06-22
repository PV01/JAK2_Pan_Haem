#################
# JAK2 V617F query ##

setwd("H:/R_projects/JAK2_Pan_Haem")

library(tidyverse)
library(odbc)
library(readxl)
library(Rsamtools)
library(IRanges)

query_db<-function(query){
  file_path<-"G:/COMPUTER/Shire/read_only_user/ForReports/db1.mdb"    
  conn <- dbConnect(drv = odbc(), .connection_string = paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=",file_path,";"))
  res<-dbGetQuery(conn, query, stringsAsFactors = FALSE)
  dbDisconnect(conn)
  return(res)
}

get_JAK2_values<-function(){
  
  q_JAK2_values<-paste("SELECT DNALAB_REPORT.LABNO, DNALAB_REPORT.INDICATION, DNALAB_REPORT.REASON, DNALAB_REPORT.SIGN_BY, DNA_Worksheet_det_result.WORKSHEET, 
                               DNA_Worksheet_det_result.LANE, DNA_Worksheet_det_result.TEST, DNA_Worksheet_det_result.RESULT, DNA_Worksheet_det_result.VALUE1, DNA_Worksheet_det_result.VALUE2, 
                               DNA_Worksheet_det_result.COMMENT
                        FROM (DNA_Worksheet_det INNER JOIN DNALAB_REPORT ON DNA_Worksheet_det.LABNO = DNALAB_REPORT.LABNO) INNER JOIN DNA_Worksheet_det_result 
                             ON (DNA_Worksheet_det.LANE = DNA_Worksheet_det_result.LANE) AND (DNA_Worksheet_det.WORKSHEET = DNA_Worksheet_det_result.WORKSHEET)
                       WHERE (((DNALAB_REPORT.LABNO) Like 'D22%' Or (DNALAB_REPORT.LABNO) Like 'D23%') AND ((DNALAB_REPORT.INDICATION)='D-MPD') AND 
                            ((DNALAB_REPORT.REASON)='Diagnosis') AND ((DNA_Worksheet_det_result.TEST)='JAK2_V617F ddPCR'));")
  
  JAK2_values<-query_db(q_JAK2_values) %>%
               filter(!grepl("No|NO", RESULT)) %>%
               mutate(RESULT=round(as.numeric(RESULT),2),
                      VALUE1=as.numeric(VALUE1),
                      VALUE2=as.numeric(VALUE2),
                      COMMENT=as.numeric(COMMENT)) %>%
              distinct(LABNO, .keep_all = TRUE)
  
  
  
  return(JAK2_values)
}


#0.25, 0.5, 1.0, 2.0 and 3.0.

JAK2_values<-get_JAK2_values()


JAK2_values<-JAK2_values %>%
                       filter((RESULT>0.22 & RESULT<0.28)|(RESULT>0.85 & RESULT<1.1)|(RESULT>1.9 & RESULT<2.1)|(RESULT>2.8 & RESULT<3.2)) %>%
                       arrange(RESULT) %>%
                       select(-SIGN_BY)

write.csv(JAK2_values,"JAK2_values.cSV")


LABNOs_selected<-JAK2_values %>%
                  pull(LABNO)



get_NGS<-function(LABNOs_selected){
  
  q_NGS<-paste("SELECT DNALAB_REPORT.LABNO, DNALAB_REPORT.INDICATION, DNALAB_REPORT.REASON, DNALAB_REPORT.SIGN_BY, DNA_Worksheet_det.WORKSHEET, 
                               DNA_Worksheet_det.LANE, DNA_Worksheet_det_result.TEST, DNA_Worksheet_det_result.RESULT, DNA_Worksheet_det_result.VALUE1, DNA_Worksheet_det_result.VALUE2, 
                               DNA_Worksheet_det_result.COMMENT
                        FROM (DNALAB_REPORT INNER JOIN DNA_Worksheet_det ON DNALAB_REPORT.LABNO = DNA_Worksheet_det.LABNO) INNER JOIN DNA_Worksheet_det_result 
                               ON (DNA_Worksheet_det.LANE = DNA_Worksheet_det_result.LANE) AND (DNA_Worksheet_det.WORKSHEET = DNA_Worksheet_det_result.WORKSHEET)
                        WHERE (((DNALAB_REPORT.LABNO) In ('",paste0(LABNOs_selected, collapse = "','"),"')) AND ((DNALAB_REPORT.INDICATION)='D-MPD') AND ((DNALAB_REPORT.REASON) Like 'NGS%') 
                                 AND ((DNA_Worksheet_det_result.TEST)='Seq NGS TSMP v3'));")
  
  NGS<-query_db(q_NGS) 
  

  return(NGS)
}

JAK2_NGS<-get_NGS(LABNOs_selected) %>%
                  select(-SIGN_BY)


JAK2_values_NGS<-JAK2_values %>%
               left_join(JAK2_NGS, by=c("LABNO"="LABNO")) %>%
               filter(!is.na(INDICATION.y))

write.csv(JAK2_values_NGS,"JAK2_values_NGS.cSV")

#################################################

NGS_file_path<-function(labnr){
      query_open<-paste("SELECT DNALAB_TEST.LABNO, DNALAB_TEST.INDICATION, DNALAB_TEST.REASON, DNALAB_TEST.TEST, DNALAB_TEST.WORKSHEET, DNALAB_TEST.LANE, DNA_Worksheet_det_result.RESULT
                     FROM DNALAB_TEST LEFT JOIN DNA_Worksheet_det_result ON (DNALAB_TEST.TEST = DNA_Worksheet_det_result.TEST) AND (DNALAB_TEST.LANE = DNA_Worksheet_det_result.LANE) 
                          AND (DNALAB_TEST.WORKSHEET = DNA_Worksheet_det_result.WORKSHEET)
               WHERE (((DNALAB_TEST.LABNO)='", paste(labnr, collapse = "', '"), "') AND ((DNALAB_TEST.REASON) Like 'NGS%') AND ((DNALAB_TEST.TEST)='Seq NGS TSMP v3'));", sep="")
        
     NGS_test<- query_db(query_open)
  
   if (nrow(NGS_test)==0|is.na(NGS_test$RESULT)){
      pat_file<-"No NGS for this patient!!!"
   } else {
     NGS_path<-"G:/DNA/RESULTS/NGS/Illumina MiSeq/HO NGS/TSMP Results"
     NGS_wks_path<-file.path(NGS_path,list.files(path=NGS_path,(glob2rx(paste0("^",substr(NGS_test$WORKSHEET,1,2),"*")))))
     NGS_test_path<-file.path(NGS_wks_path,list.files(path=NGS_wks_path,pattern=glob2rx(paste0("^",NGS_test$WORKSHEET,"*"))))
     NGS_file_pattern<-paste0("^", NGS_test$WORKSHEET,"-",NGS_test$LANE,"-",gsub("[.]","-",labnr),"*results.xlsx$")
     pat_file<-(file.path(NGS_test_path,list.files(path=NGS_test_path, pattern = glob2rx(NGS_file_pattern))))[1]
     
     igv_path<-"S:/MiSeq_data/GRCh38/TSMP_Flex/pipeline_output/"
     igv_path<-list.files(igv_path, full.names=TRUE)
     igv_path<-igv_path[grepl(substr(NGS_test$WORKSHEET,1,1), substr(igv_path, 48,60))]
     igv_path<-file.path(igv_path,NGS_test$WORKSHEET) 
     igv_path<-file.path(igv_path,list.files(igv_path,  pattern = glob2rx(paste0("*",NGS_test$WORKSHEET,"*"))))
     
     bam_path<-file.path(igv_path,paste0("alignments_TSMP_",NGS_test$WORKSHEET))
     bam_pattern<-paste0("^", NGS_test$WORKSHEET,"-",NGS_test$LANE,"-",gsub("[.]","-",labnr),"*sorted.bam$")
     bam_file<-file.path(bam_path,list.files(path=bam_path, pattern = glob2rx(bam_pattern)))
   }
     
  return(data.frame(LABNO=labnr, Excel_file=pat_file, BAM=bam_file))
     
}
     
exon_14_coverage<-function(labnr){
  pat_file<-NGS_file_path(labnr)[,2]
  exon_14_coverage<-(read_excel(pat_file, sheet= "Coverage-exon") %>%
                                 filter(Gene=="JAK2", Exon=="14"))[-c(1,2,4,22, 23)] %>%
                    mutate(LABNO=labnr)
}

NGS_coverage<-data.frame()
for (i in 1:nrow(JAK2_values_NGS)) {
  NGS_coverage<-rbind(NGS_coverage, exon_14_coverage(JAK2_values_NGS$LABNO[i]))
}
  
JAK2_values_NGS_annot<-JAK2_values_NGS %>%
                                left_join(NGS_coverage, by=c("LABNO"="LABNO"))
                                
write.csv(JAK2_values_NGS_annot,"JAK2_values_NGS_annot.cSV")

##########################################






JAK2_depth<-function(labnr){
  chrom="chr9"
  Pos=5073770
  bam_file<-NGS_file_path(labnr)[,3]
  
  bam <- scanBam(bam_file)[[1]]
  
  ranges <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))
  ranges <- GRanges(seqnames=Rle(bam$rname), ranges=ranges, strand=Rle(bam$strand), flag=bam$flag, readid=bam$rname)
  
  
  
  
  
  ranges<-
    
    ranges[seqnames=chrom] 
  
  
  
  mean(coverage(ranges))
    
    
  ranges <- IRanges(start=Pos, width=1,names=make.names(bam$qname, unique=TRUE))
  ranges <- GRanges(seqnames=Rle("chr9", ranges=ranges, strand=Rle(bam$strand), flag=bam$flag, readid=bam$rname)

                    coverage(ranges)


