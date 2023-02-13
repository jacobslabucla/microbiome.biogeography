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