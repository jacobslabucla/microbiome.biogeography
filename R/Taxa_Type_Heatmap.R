#' Generating heatmaps to match the metabolic maps containing "GMM" - gut metabolic modules
#' 
#' 
#' @author Julianne C. Yang, Jonathan P. Jacobs
#' @param dataframe dataframe containing Maaslin2 results across all six LMEM 
#' @param path_to_module_key filepath to the GMM key with annotations for each GMM
#' @param Y hierachy_annotations, feature_annotations, metabolic_map, or Map
#' @param ystring choose between "hierachy_annotations", "feature_annotations", "metabolic_map", or "Map"
#' @param titlestring title of plot, e.g. "Mucosal"
#' @param colorvector vector of colors 
#' @param breakvector vector of breaks 
#' @return a ggplot2 object encoding a heatmap
#' @export 


generate_taxa_heat_map_by_type<- function(dataframe, titlestring, colorvector, breakvector){
  #site_heatmap<-as.data.frame(heatmap_final)
  
  
  site_heatmap<-as.data.frame(dataframe)
  data<- site_heatmap
  
  
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
  
  data$Site<-factor(data$Site, levels=c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))
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
  data$Site = revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
  data$Site = factor(data$Site, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  print(summary(data$coef_d))

  data$Phylum <- gsub(".*p__", "", data$feature)
  y = tapply(data$coef_d, data$Phylum, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
  y = sort(y, FALSE)   #switch to TRUE to reverse direction
  data$Phylum= factor(as.character(data$Phylum), levels = names(y))
  
  g1 <- ggplot(data, aes(x = Site, y=Phylum)) + 
      geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
      geom_text(aes(label=asterisk)) +
      scale_fill_stepsn(breaks=bk, values = NULL, colors = cols) +
      theme_cowplot(12) +
      theme(legend.position="right",legend.justification = "center") +
      xlab("")+
      ylab("") +
      guides(fill=guide_colourbar(title="",label=TRUE,barheight = 15)) +
      ggtitle({{titlestring}}) +
      theme(plot.title = element_text(hjust = 0.5))
  
  return(g1)
}
