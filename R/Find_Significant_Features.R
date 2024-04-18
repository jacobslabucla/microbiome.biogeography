#' Generate target vector containing all significant features in all DC vs all comparisons
#' 
#' 
#' @author Julianne C. Yang, Jonathan P. Jacobs
#' @param filepath_to_significant_results_tsv a filepath to significant_results.tsv output from Maaslin2
#' @return a character vector containing all concordant and significant features across six sites
#' @export 


find_concordant_features_across_sites <- function(filepath_to_significant_results_tsv) {
  luminal<-read.table(filepath_to_significant_results_tsv, header=TRUE)
  duodenum_significant<-filter(luminal, metadata=="Site" & value=="Duodenum" &qval<0.05)
  a<-duodenum_significant$feature
  jejunum_significant<-filter(luminal, metadata=="Site" & value=="Jejunum" &qval<0.05)
  b<-jejunum_significant$feature
  ileum_significant<-filter(luminal, metadata=="Site" & value=="Ileum" &qval<0.05)
  c<-ileum_significant$feature
  cecum_significant<-filter(luminal, metadata=="Site" & value=="Cecum" &qval<0.05)
  d<-cecum_significant$feature  
  pc_significant<-filter(luminal, metadata=="Site" & value=="Proximal_Colon" &qval<0.05)
  e<-pc_significant$feature  
  DC_significant<-filter(luminal, metadata=="Site" & value=="Distal_Colon" &qval<0.05)
  f<-DC_significant$feature  
  joinab<- union(a,b)
  joincd<- union(c,d)
  joinef<- union(e,f)
  joinabcd <- union(joinab,joincd)
  target<-union(joinabcd,joinef)
  return(target)
}

#' Find concordant features across different sites
#' 
#' 
#' @author Julianne C. Yang, Jonathan P. Jacobs
#' @param target vector generated from find_features_union_for_type_heatmap
#' @return a character vector containing all concordant and significant features across six sites
#' @export 


#duo_filepath <- "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv"
#jej_filepath <- "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv"
#ile_filepath <- "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv"
#cec_filepath <- "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv"
#pc_filepath <- "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Proximal_Colon-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv"
#dc_filepath <-"Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Distal_Colon-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv"

find_features_union_for_type_heatmap <- function(duo_filepath, 
                                                 jej_filepath,
                                                 ile_filepath,
                                                 cec_filepath,
                                                 pc_filepath,
                                                 dc_filepath) {
  duodenum<-read.delim(here::here(
    path=duo_filepath))
  duodenum_significant<-filter(duodenum, metadata=="Type" & value=="Mucosal" &qval<0.05)
  a<-duodenum_significant$feature
  jejunum<-read.delim(here(
    path=jej_filepath))
  jejunum_significant<-filter(jejunum, metadata=="Type" & value=="Mucosal" &qval<0.05)
  b<-jejunum_significant$feature
  ileum<-read.delim(here(
    path=ile_filepath))
  ileum_significant<-filter(jejunum, metadata=="Type" & value=="Mucosal" &qval<0.05)
  c<-ileum_significant$feature
  cecum<-read.delim(here(
    path=cec_filepath))
  cecum_significant<-filter(cecum, metadata=="Type" & value=="Mucosal" &qval<0.05)
  d<-cecum_significant$feature  
  pc<-read.delim(here(
    path=pc_filepath))
  pc_significant<-filter(pc, metadata=="Type" & value=="Mucosal" &qval<0.05)
  e<-pc_significant$feature  
  DC<-read.delim(here(
    path=dc_filepath))
  DC_significant<-filter(DC, metadata=="Type" & value=="Mucosal" &qval<0.05)
  f<-DC_significant$feature  
  joinab<- union(a,b)
  joincd<- union(c,d)
  joinef<- union(e,f)
  joinabcd <- union(joinab,joincd)
  target<-union(joinabcd,joinef)
  
  return(target)
  
}

#' Query a vector of features against all entries within multiple dataframes
#' 
#' 
#' @author Julianne C. Yang, Jonathan P. Jacobs
#' @param target vector generated from find_features_union_for_type_heatmap
#' @return a character vector containing all concordant and significant features across six sites
#' @export 

query_type_features_union <- function(target, duo_filepath,
                                      jej_filepath,
                                      ile_filepath,
                                      cec_filepath,
                                      pc_filepath,
                                      dc_filepath) {
  duodenum<-read.table(duo_filepath, header=TRUE)
  duodenum_all<-filter(duodenum, metadata=="Type" & value=="Mucosal")
  duodenum_all<-duodenum_all[match(target,duodenum_all$feature),]
  duodenum_all$Site<- "Duodenum"
  jejunum<-read.table(jej_filepath, header=TRUE)
  jejunum_all<-filter(jejunum, metadata=="Type" & value=="Mucosal")
  jejunum_all<-jejunum_all[match(target,jejunum_all$feature),]
  jejunum_all$Site<- "Jejunum"
  ileum<-read.table(ile_filepath, header=TRUE)
  ileum_all<-filter(ileum, metadata=="Type" & value=="Mucosal")
  ileum_all<-ileum_all[match(target,ileum_all$feature),]
  ileum_all$Site<- "Ileum"
  cecum<-read.table(cec_filepath, header=TRUE)
  cecum_all<-filter(cecum, metadata=="Type" & value=="Mucosal")
  cecum_all<-cecum_all[match(target,cecum_all$feature),]
  cecum_all$Site<- "Cecum"
  pc<-read.table(pc_filepath, header=TRUE)
  pc_all<-filter(pc, metadata=="Type" & value=="Mucosal")
  pc_all<-pc_all[match(target,pc_all$feature),]
  pc_all$Site<- "Proximal_Colon"
  DC<-read.table(dc_filepath, header=TRUE)
  DC_all<-filter(DC, metadata=="Type" & value=="Mucosal")
  DC_all<-DC_all[match(target,DC_all$feature),]
  DC_all$Site<- "Distal_Colon"
  
  duojej<-rbind(duodenum_all,jejunum_all)
  ilecec<-rbind(ileum_all, cecum_all)
  pcdc<-rbind(pc_all,DC_all)
  duojejilecec<-rbind(duojej,ilecec)
  duojejilececpcdc<-rbind(duojejilecec,pcdc)
  return(duojejilececpcdc)
  
}
