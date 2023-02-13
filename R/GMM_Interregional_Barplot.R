#' Generating heatmaps to match the metabolic maps containing "GMM" - gut metabolic modules
#' 
#' 
#' @author Julianne C. Yang, Jonathan P. Jacobs
#' @param path_to_all_results_tsv filepath to Maaslin2/all_results.tsv
#' @param path_to_module_key filepath to the GMM key with annotations for each GMM
#' @param ystring choose between "hierachy_annotations", "feature_annotations", "metabolic_map", or "Map"
#' @param titlestring title of plot, e.g. "Mucosal"
#' @param colorvector character vector of two colors
#' @return a ggplot2 object encoding a heatmap
#' @export 


generate_interregional_GMM_barplot <- function(path_to_all_results_tsv, path_to_Module_Key, ystring,titlestring, colorvector){
  ##### test function inputs 
  #annotation <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv", header=TRUE)
  #luminal <- read.delim("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv")
  
  #cols <- viridis::viridis(2)
  #####
  
  luminal<-read.table(path_to_all_results_tsv, header=TRUE)
  luminal <- luminal %>% filter(metadata=="Site_General" & qval<0.05)
  
  annotation <- read.csv(path_to_Module_Key, header=TRUE)
  
  data<- (merge(luminal, annotation, by = 'feature'))
  cols=c(colorvector)

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
    
    res_plot <- data %>% select(c("coef", "qval","hierachy_annotations"))
    res_plot <- unique(res_plot)
    res_plot <- res_plot %>%
      mutate(site = ifelse(coef< 0, "Colon", "SI"))
    
      y = tapply(res_plot$coef, res_plot$hierachy_annotations, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
      y = sort(y, FALSE)   #switch to TRUE to reverse direction
      res_plot$hierachy_annotations= factor(as.character(res_plot$hierachy_annotations), levels = names(y))
     
    g1<- res_plot %>%
      arrange(coef) %>%
      filter(qval < 0.05, abs(coef) > 0) %>%
      ggplot2::ggplot(aes(coef, hierachy_annotations, fill = site)) +
        geom_bar(stat = "identity") +
        cowplot::theme_cowplot(12) +
        theme(axis.text.y = element_text(face = "italic")) +
        scale_fill_manual(values = cols) +
        labs(x = "Effect size (SI/Colon)",
           y = "",
           fill = "") +
      theme(legend.position = "none")+
      ggtitle({{titlestring}}) +
      theme(plot.title = element_text(hjust = 0.5))
    
  }
  else if(ystring=="metabolic_map"){
    res_plot <- data %>% select(c("coef", "qval","metabolic_map"))
    res_plot <- unique(res_plot)
    res_plot <- res_plot %>%
      mutate(site = ifelse(coef< 0, "Colon", "SI"))
    
    y = tapply(res_plot$coef, res_plot$metabolic_map, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    res_plot$metabolic_map= factor(as.character(res_plot$metabolic_map), levels = names(y))
    
    g1<- res_plot %>%
      arrange(coef) %>%
      filter(qval < 0.05, abs(coef) > 0) %>%
      ggplot2::ggplot(aes(coef, metabolic_map, fill = site)) +
      geom_bar(stat = "identity") +
      cowplot::theme_cowplot(12) +
      theme(axis.text.y = element_text(face = "italic")) +
      scale_fill_manual(values = cols) +
      labs(x = "Effect size (SI/Colon)",
           y = "",
           fill = "") +
      theme(legend.position = "none")+
      ggtitle({{titlestring}}) +
      theme(plot.title = element_text(hjust = 0.5))
    
  }
  else if(ystring=="feature_annotations"){
    res_plot <- data %>% select(c("coef", "qval","feature_annotations"))
    res_plot <- unique(res_plot)
    res_plot <- res_plot %>%
      mutate(site = ifelse(coef< 0, "Colon", "SI"))
    
    y = tapply(res_plot$coef, res_plot$feature_annotations, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    res_plot$feature_annotations= factor(as.character(res_plot$feature_annotations), levels = names(y))
    
    g1<- res_plot %>%
      arrange(coef) %>%
      filter(qval < 0.05, abs(coef) > 0) %>%
      ggplot2::ggplot(aes(coef, feature_annotations, fill = site)) +
      geom_bar(stat = "identity") +
      cowplot::theme_cowplot(12) +
      theme(axis.text.y = element_text(face = "italic")) +
      scale_fill_manual(values = cols) +
      labs(x = "Effect size (SI/Colon)",
           y = "",
           fill = "") +
      theme(legend.position = "none")+
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
    
    res_plot <- ggplotdata %>% select(c("coef_mean", "qval_d","Map"))
    res_plot <- unique(res_plot)
    res_plot <- res_plot %>%
      mutate(site = ifelse(coef_mean< 0, "Colon", "SI"))
    
    y = tapply(res_plot$coef_mean, res_plot$Map, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    res_plot$Map= factor(as.character(res_plot$Map), levels = names(y))
    
    g1<- res_plot %>%
      arrange(coef_mean) %>%
      filter(qval_d < 0.05, abs(coef_mean) > 0) %>%
      ggplot2::ggplot(aes(coef_mean, Map, fill = site)) +
      geom_bar(stat = "identity") +
      cowplot::theme_cowplot(12) +
      theme(axis.text.y = element_text(face = "italic")) +
      scale_fill_manual(values = cols) +
      labs(x = "Effect size (SI/Colon)",
           y = "",
           fill = "") +
      theme(legend.position = "none")+
      ggtitle({{titlestring}}) +
      theme(plot.title = element_text(hjust = 0.5))
    
    
  }
  
  return(g1)
}
