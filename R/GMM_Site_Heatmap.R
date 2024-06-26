#' Generating heatmaps to match the metabolic maps containing "GMM" - gut metabolic modules
#' 
#' 
#' @author Julianne C. Yang, Jonathan P. Jacobs
#' @param path_to_all_results_tsv filepath to Maaslin2/all_results.tsv
#' @param targetvector vector containing the featureIDs that you want to include in the heatmap. these featureIDs must also be in all_results.tsv
#' @param path_to_module_key filepath to the GMM key with annotations for each GMM
#' @param Y hierachy_annotations, feature_annotations, metabolic_map, or Map
#' @param ystring choose between "hierachy_annotations", "feature_annotations", "metabolic_map", or "Map"
#' @param titlestring title of plot, e.g. "Mucosal"
#' @return a ggplot2 object encoding a heatmap
#' @export 


generate_GMM_heat_map_by_site <- function(path_to_all_results_tsv, targetvector, path_to_Module_Key, Y, ystring,titlestring, colorvector,breakvector){
  #targetvector<-muctarget
  #annotation <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv", header=TRUE)
  #luminal<-read.table("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv", header=TRUE)
  
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
  
  
  #construct the heatmap using ggplot
  annotation <- read.csv(path_to_Module_Key, header=TRUE)
  
  data<- (merge(site_heatmap, annotation, by = 'feature'))
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
  data$value = plyr::revalue(data$value, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
  data$value = factor(data$value, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  print(summary(data$coef_d))
  
  #Add ammonia and carbon dioxide and remove lactaldehyde degradation
  data <- data[!grepl("#N/A",data$Map), ]
  
  if("Map2_ammonia" %in% colnames(data) | "Map3_carbon_dioxide" %in% colnames(data)){
    if("Map2_ammonia" %in% colnames(data)) {
      ammonia<- data %>% filter(Map2_ammonia=="ammonia") 
      if(nrow(ammonia)!=0) {
        ammonia$Map <- "ammonia"
      }
      data<-rbind(data, ammonia)
    }
    if("Map3_carbon_dioxide" %in% colnames(data)){
      co2 <- data %>% filter(Map3_carbon_dioxide=="carbon dioxide")
      if(nrow(co2)!=0) {
        co2$Map <- "carbon dioxide"
      }
      data <- rbind(data, co2)
    }
  }
  
  #Generate Y values (what you want to show on Y axis)
  data$feature_annotations<-paste(data$feature,data$annotation,sep=" : ")
  data$hierachy_annotations<-paste(data$Hierarchy_L2,data$annotation,sep=" : ")
  data$metabolic_map<-paste(data$Map,data$annotation,sep=" : ")
  
  if(ystring=="hierachy_annotations"){
    y = tapply(data$coef, data$hierachy_annotations, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    data$hierachy_annotations= factor(as.character(data$hierachy_annotations), levels = names(y))
    ggplotdata<-data
    
    g1 <- ggplot(ggplotdata, aes(x = value, y={{Y}})) + 
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
  else if(ystring=="metabolic_map"){
    y = tapply(data$coef_d, data$metabolic_map, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    data$metabolic_map= factor(as.character(data$metabolic_map), levels = names(y))
    ggplotdata<-data
    
    g1 <- ggplot(ggplotdata, aes(x = value, y={{Y}})) + 
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
  else if(ystring=="feature_annotations"){
    y = tapply(data$coef_d, data$feature_annotations, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    data$feature_annotations= factor(as.character(data$feature_annotations), levels = names(y))
    ggplotdata<-data
  
    g1 <- ggplot(ggplotdata, aes(x = value, y={{Y}})) + 
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
  
  else if(ystring=="Map"){
    
    ggplotdata<-data
    
    ggplotdata1<- ggplotdata %>%
      group_by(value,Map) %>%
      summarise_at(vars(qval), list(qval_d = min))
    ggplotdata1$Site_Map <- paste(ggplotdata1$value, ggplotdata1$Map)
    
    ggplotdata2<- ggplotdata %>%
      group_by(value,Map) %>%
      summarise_at(vars(coef), list(coef_mean = median))
    ggplotdata2$Site_Map <- paste(ggplotdata2$value, ggplotdata2$Map)
    
    ggplotdata <- merge(ggplotdata1, ggplotdata2, by= "Site_Map")
    ggplotdata$Map <- ggplotdata$Map.x
    ggplotdata$value <- ggplotdata$value.x
    
    qval<-ggplotdata$qval_d
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
    ggplotdata$asterisk<-asterisk
    
    
    ggplotdata <- ggplotdata %>% mutate(coef_d_mean= ifelse(coef_mean>1.5, 2, coef_mean))
    ggplotdata$coef_d_mean[ggplotdata$coef_d_mean < (-1.5)] <- (-2)
    if(max(ggplotdata$coef_d_mean) >1.0 & max(ggplotdata$coef_d_mean)<1.5){
    ggplotdata$coef_d_mean[ggplotdata$coef_d_mean >= 1.0 & ggplotdata$coef_d_mean<(1.5)] <- (1.5)
    }
    if(max(ggplotdata$coef_d_mean) >0 & max(ggplotdata$coef_d_mean)<1.0){
      ggplotdata$coef_d_mean[ggplotdata$coef_d_mean >= 0.5 & ggplotdata$coef_d_mean<1.0] <- (1.0)
    }
    if(max(ggplotdata$coef_d_mean) >0 & max(ggplotdata$coef_d_mean)<0.5){
      ggplotdata$coef_d_mean[ggplotdata$coef_d_mean >= 0 & ggplotdata$coef_d_mean<0.5] <- (0.5)
    }
    if(min(ggplotdata$coef_d_mean) < (-1.0) & min(ggplotdata$coef_d_mean)> (-1.5)){
      ggplotdata$coef_d_mean[ggplotdata$coef_d_mean < (-1.0) & ggplotdata$coef_d_mean>=(-1.5)] <- (-1.5)
    }
    if(min(ggplotdata$coef_d_mean) < (-0.5) & min(ggplotdata$coef_d_mean)>(-1.0)){
      ggplotdata$coef_d_mean[ggplotdata$coef_d_mean < (-0.5) & ggplotdata$coef_d_mean>=(-1.0)] <- (-1.0)
    }
    if(min(ggplotdata$coef_d_mean) < (0) & min(ggplotdata$coef_d_mean)>(-0.5)){
      ggplotdata$coef_d_mean[ggplotdata$coef_d_mean < (0) & ggplotdata$coef_d_mean>=(-0.5)] <- (0)
    }
   
  
    
    
    y = tapply(ggplotdata$coef_d_mean, ggplotdata$Map, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    ggplotdata$Map= factor(as.character(ggplotdata$Map), levels = names(y))
  
    print(summary(ggplotdata$coef_d_mean))
    
    g1 <- ggplot(ggplotdata, aes(x = value, y={{Y}})) + 
      geom_tile(aes(fill = coef_d_mean),colour="white",size=0.25) +
      geom_text(aes(label=asterisk)) +
      #scale_fill_stepsn(colors = cols) +
      scale_fill_stepsn(breaks=bk, values = NULL, colors = cols) +
      cowplot::theme_cowplot(12) +
      theme(legend.position="top",legend.justification = "center") +
      xlab("")+
      ylab("") +
      guides(fill=guide_colourbar(title="",label=TRUE,barwidth = 15)) +
      ggtitle({{titlestring}}) +
      theme(plot.title = element_text(hjust = 0.5))
    
    
  }
  
  return(g1)
}
