#' Generate taxa summary plots by genera showing only the top 0.1% most abundant 
#'
#' 
#' Plot does not have legend - legends were graphed separately 
#' 
#' 
#' @author Julianne Yang
#' @param path_to_RDS filepath to where your L6 RDS is stored
#' @param titlestring title of plot
#' @param greppattern ".*s__" for species 
#' @param fillvector a named vector of colors 
#' @param graphby a string, can pass "Site" or "Type"
#' @return a ggplot2 object which you can then further customize by adding ggplot2 functions
#' @export
#' @examples
#' 
#' generate_L6_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Luminal_level-6.RDS", 
#' "Luminal ( > 0.1% Relative Abundance)",
#'  ".*g__",
#'   assign_cols,
#'    "Site")
#' 
#' 

generate_L6_taxa_plots <- function(path_to_RDS, titlestring,greppattern, fillvector, graphby){
  titlestring<-c(titlestring)
  L2_lum<-readRDS(path_to_RDS)
  L2_lum<- as.matrix(L2_lum)
  L2_lum<-funrar::make_relative(L2_lum)
  L2_lum<-as.data.frame(t(L2_lum))
  toptaxa<- rowMeans(L2_lum)
  L2_lum$averageRA <-toptaxa/6
  L2_lum <- L2_lum %>% dplyr::mutate(keeptaxa = ifelse(averageRA >0.001, row.names(L2_lum), "Other"))
  L2_lum <-select(L2_lum,-averageRA)
  
  taxa<-L2_lum$keeptaxa
  L2_lum <- select(L2_lum,-keeptaxa)
  L2_lum <- as.matrix(sapply(L2_lum,as.numeric))
  L2_lum <- as.data.frame(prop.table(L2_lum,2))
  taxa<-gsub(greppattern,"",taxa )
  
  L2_lum$Taxa <-taxa
  L2_lum<- tidyr::pivot_longer(L2_lum, -c(Taxa), values_to ="Value", names_to ="Site")
  L2_lum$Value <- L2_lum$Value * 100
  
  if({{graphby}} == "Site"){
    L2_lum$Site = plyr::revalue(L2_lum$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
    L2_lum$Site = factor(L2_lum$Site, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  }
  else if ({{graphby}}=="Type"){ 
    L2_lum$Site = plyr::revalue(L2_lum$Site, c("Luminal"="Lum", "Mucosal" = "Muc"))
  }
  
  cols <- fillvector
  ggplot2::ggplot(data=L2_lum, aes(x=Site, y=Value, fill=Taxa)) +
    geom_bar(stat="identity")+
    #scale_fill_paletteer_d(palette="colorBlindness::SteppedSequential5Steps") +
    #scale_fill_paletteer_d(palette="dutchmasters::milkmaid") +
    #scale_fill_paletteer_d("tvthemes::rickAndMorty")+
    #scale_fill_paletteer_d("ggsci::category20_d3")+
    scale_fill_manual(values = cols)+
    theme(legend.position = "none")+
    theme_cowplot(12) +
    ylab("% Relative Abundance") +
    xlab("")+
    labs(fill="") +
    ggtitle(titlestring) +
    theme(legend.position="top") +
    theme(plot.title = element_text(hjust = 0.5))+
    #guides(fill=guide_legend(nrow=8, byrow=TRUE)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
}

generate_L2_taxa_plots <- function(input_data, titlestring,greppattern, graphby){
  titlestring<-c(titlestring)
  L2_lum<-read.csv(input_data)
  L2_lum <- as.data.frame(t(L2_lum))
  colnames(L2_lum)<- L2_lum[1,]
  L2_lum <- L2_lum[-1,]
  taxa<-row.names(L2_lum)
  L2_lum <- as.matrix(sapply(L2_lum,as.numeric))
  L2_lum <- as.data.frame(prop.table(L2_lum,2))
  taxa<-gsub(greppattern,"",taxa )
  
  L2_lum$Taxa <-taxa
  L2_lum<- tidyr::pivot_longer(L2_lum, -c(Taxa), values_to ="Value", names_to ="Site")
  if({{graphby}} == "Site"){
    L2_lum$Site = plyr::revalue(L2_lum$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
    L2_lum$Site = factor(L2_lum$Site, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  }
  else if ({{graphby}}=="Type"){ 
    L2_lum$Site = plyr::revalue(L2_lum$Site, c("Luminal"="Lum", "Mucosal" = "Muc"))
  }
  L2_lum$Value <- L2_lum$Value * 100
  
  
  ggplot2::ggplot(data=L2_lum, aes(x=Site, y=Value, fill=Taxa)) +
    geom_bar(stat="identity")+
    scale_fill_viridis_d(option="B")+
    theme(legend.position = "right")+
    theme_cowplot(12) +
    ylab("% Relative Abundance") +
    xlab("")+
    labs(fill="")+
    ggtitle({{titlestring}}) +
    theme(legend.position="top") +
    theme(plot.title = element_text(hjust = 0.5))+
    guides(fill=guide_legend(nrow=4, byrow=TRUE))
  
  
}