### The birth of package "Microbiome- Biogeography" ---

generate_GMM_heat_map_by_site <- function(path_to_all_results_tsv, targetvector, path_to_Module_Key, Y, ystring,titlestring){
  
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
  data$feature_annotations<-paste(data$feature,data$annotation,sep=" : ")
  data$hierachy_annotations<-paste(data$Hierarchy_L2,data$annotation,sep=" : ")
  data$metabolic_map<-paste(data$Map,data$annotation,sep=" : ")
  
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
  data <- data %>% mutate(coef_d= ifelse(coef>2, 2, coef))
  data$coef_d[data$coef_d < (-2)] <- (-2)
  data$value = revalue(data$value, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
  data$value = factor(data$value, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  
  if(ystring=="hierachy_annotations"){
    y = tapply(data$coef, data$hierachy_annotations, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    data$hierachy_annotations= factor(as.character(data$hierachy_annotations), levels = names(y))
  }
  else if(ystring=="metabolic_map"){
    y = tapply(data$coef_d, data$metabolic_map, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    data$metabolic_map= factor(as.character(data$metabolic_map), levels = names(y))
  }
  else if(ystring=="feature_annotations"){
    y = tapply(data$coef_d, data$feature_annotations, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    data$feature_annotations= factor(as.character(data$feature_annotations), levels = names(y))
  }
  
  ggplotdata<-data
  cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
  bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
  
  g1 <- ggplot(ggplotdata, aes(x = value, y={{Y}})) + 
    geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
    geom_text(aes(label=asterisk)) +
    scale_fill_stepsn(breaks=bk, values = NULL, colors = cols) +
    theme_cowplot(12) +
    theme(legend.position="top",legend.justification = "center") +
    xlab("")+
    ylab("") +
    guides(fill=guide_colourbar(title="",label=TRUE,barwidth = 15)) +
    ggtitle({{titlestring}}) +
    theme(plot.title = element_text(hjust = 0.5))
  
  g1
}