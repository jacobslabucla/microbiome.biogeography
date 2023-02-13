#' Generating heatmaps for "GBM" - gut brain modules
#' 
#' 
#' @author Julianne C. Yang, Jonathan P. Jacobs
#' @param path_to_all_results_tsv filepath to Maaslin2/all_results.tsv
#' @param targetvector vector containing the featureIDs that you want to include in the heatmap. these featureIDs must also be in all_results.tsv
#' @param path_to_module_key filepath to the GMM key with annotations for each GMM
#' @param titlestring title of plot, e.g. "Mucosal"
#' @return a ggplot2 object encoding a heatmap
#' @export 

find_features_union_for_type_heatmap <- function(duo_filepath, 
                                                 jej_filepath,
                                                 ile_filepath,
                                                 cec_filepath,
                                                 pc_filepath,
                                                 dc_filepath) {
duodenum<-readr::read_delim(here::here(
  path=duo_filepath)
  ,col_types = list(.default = readr::col_guess()))
duodenum_significant<-filter(duodenum, metadata=="Type" & value=="Mucosal" &qval<0.05)
a<-duodenum_significant$feature
jejunum<-readr::read_delim(here::here(
  path=jej_filepath)
  ,col_types = list(.default = readr::col_guess()))
jejunum_significant<-filter(jejunum, metadata=="Type" & value=="Mucosal" &qval<0.05)
b<-jejunum_significant$feature
ileum<-readr::read_delim(here::here(
  path=ile_filepath)
  ,col_types = list(.default = readr::col_guess()))
ileum_significant<-filter(jejunum, metadata=="Type" & value=="Mucosal" &qval<0.05)
c<-ileum_significant$feature
cecum<-readr::read_delim(here::here(
  path=cec_filepath)
  ,col_types = list(.default = readr::col_guess()))
cecum_significant<-filter(cecum, metadata=="Type" & value=="Mucosal" &qval<0.05)
d<-cecum_significant$feature  
pc<-readr::read_delim(here::here(
  path=pc_filepath),
  col_types = list(.default = readr::col_guess()))
pc_significant<-filter(pc, metadata=="Type" & value=="Mucosal" &qval<0.05)
e<-pc_significant$feature  
DC<-readr::read_delim(here::here(
  path=dc_filepath),
  col_types = list(.default = readr::col_guess()))
DC_significant<-filter(DC, metadata=="Type" & value=="Mucosal" &qval<0.05)
f<-DC_significant$feature  
joinab<- union(a,b)
joincd<- union(c,d)
joinef<- union(e,f)
joinabcd <- union(joinab,joincd)
target<-union(joinabcd,joinef)

return(target)

}

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

generate_taxa_heat_map_by_site <- function(path_to_all_results_tsv, targetvector, titlestring, colorvector,breakvector){
  #targetvector<-lumtarget
  #luminal<-read.table("Maaslin2_L2/UCLA_O_SPF/L2-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv", header=TRUE)
  #cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
  #bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
  
  luminal<-read.table(path_to_all_results_tsv, header=TRUE)
  
  target <- targetvector
  luminal_all<-filter(luminal, metadata=="Site")
  data<-luminal_all[luminal_all$feature %in% target, ]
  
  # make an empty dataframe to store the reference variable and assign x, a string vector, to y as its column names:
  y <- data.frame(matrix(NA,nrow=length(target),ncol=9))
  x <- c(colnames(data))
  colnames(y) <- x
  y$feature<-target
  y$coef <- 0
  y$value <- "Distal_Colon"
  y$metadata <-"Site"
  y$qval<-100
  
  site_heatmap<-rbind(data,y)
  site_heatmap <- site_heatmap %>% select(c("metadata","feature","value","coef","qval"))
  
  #construct the heatmap using ggplot
  site_heatmap$Phylum <- gsub(".*p__","",site_heatmap$feature)
  data<- site_heatmap
  
  cols=c(colorvector)
  bk =c(breakvector)
  
  #write.csv(data,"Luminal-DCvsall-ggplot-Heatmap.csv")
  qval<-data$qval
  asterisk<-c("")
  for (item in qval){
    if (item < 0.05){
      asterisk<-c(asterisk,"*")
    }
    else {
      asterisk<-c(asterisk,"")
    }
  }
  asterisk<-asterisk[-1]
  data$asterisk<-asterisk
  data$value<-factor(data$value, levels=c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))
  data <- data %>% mutate(coef_d= ifelse(coef>1.5, 2, coef))
  data$coef_d[data$coef_d < (-1.5)] <- (-2)
  if(max(data$coef_d) >1.0 & max(data$coef_d)<1.5){
    data$coef_d[data$coef_d >= 1.0 & data$coef_d<(1.5)] <- (1.5)
  }
  if(max(data$coef_d) >0 & max(data$coef_d)<1.0){
    data$coef_d[data$coef_d >= 0.5 & data$coef_d<1.0] <- (1.0)
  }
  if(max(data$coef_d) >0 & max(data$coef_d)<0.5){
    data$coef_d[data$coef_d >= 0 & data$coef_d<0.5] <- (0.5)
  }
  if(min(data$coef_d) < (-1.0) & min(data$coef_d)> (-1.5)){
    data$coef_d[data$coef_d < (-1.0) & data$coef_d>=(-1.5)] <- (-1.5)
  }
  if(min(data$coef_d) < (-0.5) & min(data$coef_d)>(-1.0)){
    data$coef_d[data$coef_d < (-0.5) & data$coef_d>=(-1.0)] <- (-1.0)
  }
  if(min(data$coef_d) < (0) & min(data$coef_d)>(-0.5)){
    data$coef_d[data$coef_d < (0) & data$coef_d>=(-0.5)] <- (0)
  }
  data$value = revalue(data$value, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
  data$value = factor(data$value, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  print(summary(data$coef_d))
  
  data$annotation <- data$Phylum
  y = tapply(data$coef_d, data$annotation, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
  y = sort(y, FALSE)   #switch to TRUE to reverse direction
  data$annotation= factor(as.character(data$annotation), levels = names(y))
  ggplotdata<-data
  
  ggplot(ggplotdata, aes(x = value, y=annotation)) + 
    geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
    geom_text(aes(label=asterisk)) +
    scale_fill_stepsn(breaks=bk, values = NULL, colors = cols) +
    cowplot::theme_cowplot(12) +
    theme(legend.position="top",legend.justification = "center") +
    xlab("")+
    ylab("") +
    guides(fill=guide_colourbar(title="",label=TRUE,barwidth = 15)) +
    ggtitle({{titlestring}}) +
    theme(plot.title = element_text(hjust = 0.5))
  
}
