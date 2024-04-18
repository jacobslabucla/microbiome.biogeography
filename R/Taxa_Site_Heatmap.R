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
